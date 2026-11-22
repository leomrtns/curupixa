#!/usr/bin/env python3
"""Optional SciPy/scikit-learn comparisons; pass a shared phylogenetics library as argv[1].
The HDBSCAN oracle cuts a SciPy MST at each distinct weight, independently of Curupixa's merge-tree construction.
Requires numpy, scipy and scikit-learn in the selected Python environment; none are Curupixa dependencies.
"""
import ctypes as ct
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from sklearn.cluster import HDBSCAN
from sklearn.metrics import adjusted_rand_score


class DensityResult(ct.Structure):
    _fields_ = [("samples", ct.c_size_t), ("clusters", ct.c_size_t), ("edge_count", ct.c_size_t),
                ("labels", ct.POINTER(ct.c_size_t)), ("selected", ct.c_void_p),
                ("probabilities", ct.POINTER(ct.c_double)), ("stabilities", ct.c_void_p), ("hierarchy", ct.c_void_p)]


def threshold_oracle(matrix, min_points, minimum, method, allow_root):
    """Connected components below each MST cut define simultaneous splits, including tied core distances."""
    core = np.partition(matrix, min_points - 1, axis=1)[:, min_points - 1]
    reachability = np.maximum(matrix, np.maximum(core[:, None], core[None, :]))
    np.fill_diagonal(reachability, 0)
    tree = minimum_spanning_tree(reachability).tocoo()
    edges = list(zip(tree.row.tolist(), tree.col.tolist(), tree.data.tolist()))

    def condense(members, birth):
        original = members.copy()
        departures = {}
        children = []
        while True:
            internal = [(a, b, w) for a, b, w in edges if a in members and b in members]
            weight = max(w for a, b, w in internal)
            density = np.longdouble(1) / np.longdouble(weight) if weight else np.longdouble(np.inf)
            adjacency = {point: [] for point in members}
            for a, b, w in internal:
                if w < weight:
                    adjacency[a].append(b)
                    adjacency[b].append(a)
            unseen = members.copy()
            groups = []
            while unseen:
                group = {unseen.pop()}
                stack = list(group)
                while stack:
                    for point in adjacency[stack.pop()]:
                        if point in unseen:
                            unseen.remove(point)
                            group.add(point)
                            stack.append(point)
                groups.append(group)
            large = [group for group in groups if len(group) >= minimum]
            if len(large) == 1:
                for point in members - large[0]:
                    departures[point] = density
                members = large[0]
            else:
                departures.update({point: density for point in members})
                children = [condense(group, density) for group in large]
                break
        stability = sum((0 if value == birth else value - birth) for value in departures.values())
        return {"members": original, "departures": departures, "stability": stability, "children": children}

    root = condense(set(range(len(matrix))), np.longdouble(0))

    def select(node, is_root=False):
        selected = [cluster for child in node["children"] for cluster in select(child)]
        eligible = not is_root or allow_root
        if method == "leaf":
            return [node] if not is_root and not node["children"] else selected
        if eligible and node["stability"] > 0 and node["stability"] >= sum(c["stability"] for c in selected):
            return [node]
        return selected

    labels = np.full(len(matrix), -1, dtype=int)
    probabilities = np.zeros(len(matrix))
    for label, cluster in enumerate(select(root, True)):
        maximum = max(cluster["departures"].values())
        for point, value in cluster["departures"].items():
            if cluster is root and value < maximum:
                continue
            labels[point] = label
            probabilities[point] = 1 if np.isinf(value) or maximum == 0 else value / maximum
    return labels, probabilities


def main():
    library = ct.CDLL(str(Path(sys.argv[1]).resolve()))
    pointer, size, real = ct.c_void_p, ct.c_size_t, ct.c_double
    callback_type = ct.CFUNCTYPE(ct.c_int, pointer, size, size, pointer, size)
    library.crpx_assignment_jv.argtypes = [pointer, size, size, pointer, pointer]
    library.crpx_distance_generator_new.argtypes = [size, size, callback_type, pointer, ct.c_bool, pointer]
    library.crpx_distance_generator_free.argtypes = [pointer]
    library.crpx_distance_generator_free.restype = None
    library.crpx_hdbscan_run.argtypes = [pointer, size, size, size, ct.c_int, ct.c_bool, pointer]
    library.crpx_hdbscan_free.argtypes = [pointer]
    library.crpx_hdbscan_free.restype = None
    rng = np.random.default_rng(728)
    for _ in range(500):
        rows = int(rng.integers(1, 31))
        columns = int(rng.integers(rows, rows + 20))
        costs = np.ascontiguousarray(rng.uniform(0, 100, (rows, columns)))
        if _ % 2:
            costs = np.floor(costs)
        assignment = np.zeros(rows, dtype=np.uintp)
        total = real()
        assert library.crpx_assignment_jv(costs.ctypes.data, rows, columns, assignment.ctypes.data, ct.byref(total)) == 0
        a, b = linear_sum_assignment(costs)
        assert np.isclose(total.value, costs[a, b].sum(), atol=1e-10)
    print("500 assignments agree with SciPy", flush=True)
    checked = sklearn_exact = 0
    for iteration in range(120):
        count, dimensions = int(rng.integers(10, 65)), int(rng.integers(1, 5))
        points = rng.normal(size=(count, dimensions))
        points[count // 2:] += 4
        if iteration % 3 == 0:
            points[-1] = 20
        matrix = np.ascontiguousarray(cdist(points, points))

        @callback_type
        def distances(_data, first, second, out, components):
            ct.cast(out, ct.POINTER(real))[0] = matrix[first, second]
            return 0

        generator = pointer()
        assert library.crpx_distance_generator_new(count, 1, distances, None, False, ct.byref(generator)) == 0
        min_points = 1 if iteration % 2 == 0 else int(rng.integers(2, min(10, count)))
        minimum = int(rng.integers(2, min(10, count)))
        for method in ["eom", "leaf"]:
            for allow_root in [False, True]:
                result = ct.POINTER(DensityResult)()
                assert library.crpx_hdbscan_run(generator, 0, min_points, minimum, int(method == "leaf"),
                                                allow_root, ct.byref(result)) == 0
                labels = np.array([int(result.contents.labels[i]) if result.contents.labels[i] != size(-1).value else -1
                                   for i in range(count)])
                probabilities = np.array([result.contents.probabilities[i] for i in range(count)])
                expected, membership = threshold_oracle(matrix, min_points, minimum, method, allow_root)
                assert adjusted_rand_score(labels, expected) == 1 and np.array_equal(labels == -1, expected == -1)
                assert np.allclose(probabilities, membership, rtol=1e-10, atol=1e-10)
                oracle = HDBSCAN(min_samples=min_points, min_cluster_size=minimum, metric="precomputed", copy=True,
                                 cluster_selection_method=method, allow_single_cluster=allow_root).fit(matrix)
                exact = adjusted_rand_score(labels, oracle.labels_) == 1 and np.array_equal(labels == -1, oracle.labels_ == -1)
                exact = exact and np.allclose(probabilities, oracle.probabilities_, rtol=1e-10, atol=1e-10)
                if exact:
                    sklearn_exact += 1
                if min_points == 1:
                    assert exact, (iteration, min_points, minimum, method, allow_root, labels, oracle.labels_, probabilities, oracle.probabilities_)
                checked += 1
                library.crpx_hdbscan_free(result)
        library.crpx_distance_generator_free(generator)
    print(f"{checked} HDBSCAN configurations agree with the independent simultaneous-cut oracle", flush=True)
    print(f"{sklearn_exact}/{checked} also match scikit-learn exactly; simultaneous tied splits may change boundary labels or membership normalisation")


if __name__ == "__main__":
    main()
