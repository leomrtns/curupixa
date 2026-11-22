# Root selection from branch lengths
`rooting.h` provides midpoint, minimum-variance and minimal ancestor deviation (MAD) root selection in `curupixa_phylogenetics`. Each function reads a tree and returns an edge, a position along that edge and the value of its objective. Root selection leaves the tree unchanged. A separate call to `crpx_tree_reroot_edge()` applies the selected location.
The stored root supplies an orientation for traversal. These calculations search the entire metric tree, treating an existing binary root as a subdivision of an edge. The resulting location is an estimate of the biological root under the assumptions of the selected method. Applying that location preserves the topology and distances between tips; it does not impose a clock or make the branch lengths ultrametric.
## Input and output
Supply a structurally valid tree with at least two tips and complete, finite, non-negative lengths on every edge. The root stem is ignored. Multifurcations and unary internal vertices are accepted. A unary storage root is rejected with `CRPX_ERROR_UNSUPPORTED`, because traversing its stem would add a terminal vertex without a taxon. Remove that stem or supply a different representation before selecting a root.
A `crpx_root_location` contains `child`, `fraction` and `score`. The selected edge joins `child` to its parent in the input tree. If its length is $l$, the returned position is $x=l\,\mathtt{fraction}$ from the parent. Fractions zero and one denote existing vertices; a zero-length edge is represented with fraction zero. The location refers to the input node IDs and orientation, so calculate it again after editing the tree.
The methods minimise different quantities. Midpoint returns the maximum distance from its root to a tip, in the input length units. Minimum variance returns the population variance of root-to-tip distances, in squared input units. MAD returns a dimensionless root-mean-square deviation between zero and one. These scores describe clock-likeness according to their respective definitions; they are not posterior root probabilities.
Exactly equal calculated scores are resolved by the smallest child ID, with the smallest fraction on a constant edge. Floating-point rounding can resolve mathematically tied candidates differently after changing the representation. Several adjacent edges can represent the same optimal vertex, and a connected set of zero-length edges can represent the same metric point. The interface returns one representative rather than enumerating tied locations.
## Using a selected location
```c
#include <rooting.h>

crpx_root_location location;
crpx_status status = crpx_tree_root_mad (tree, NULL, &location);
if (status == CRPX_OK) {
  crpx_node_id root;
  status = crpx_tree_reroot_edge (tree, location.child, location.fraction,
                                 CRPX_TREE_KEEP_OLD_ROOT, &root);
}
```
`CRPX_TREE_KEEP_OLD_ROOT` preserves the old vertex and appends a new root; the old vertex can become unary. For an existing binary root in a tree with more than two tips, `CRPX_TREE_REUSE_BINARY_ROOT` instead suppresses that subdivision and reuses its node ID. Choose the policy according to how internal labels and external arrays should be handled. The existing rerooting operation moves the old root label when reusing it; reassess support annotations after rerooting.
Passing `NULL` as the calculation workspace allocates temporary storage. For repeated calculations, construct a `crpx_rooting_workspace` with `crpx_rooting_workspace_new()`, pass it to any of the three methods, and release it with `crpx_rooting_workspace_free()`. Capacity is retained between calls, including calls on unrelated trees. No tree pointer is retained. Concurrent calculations may share a read-only tree when each uses a separate workspace. Outputs remain unchanged on failure, and no worker threads or global mutable state are introduced.
## Midpoint
Let $d(u,v)$ be the sum of edge lengths on the path between two points. Find a pair of tips $a,b$ attaining the diameter
$$
D=\max_{i,j}d(i,j).
$$
The midpoint root lies halfway along their path. Equivalently, it minimises
$$
M(r)=\max_i d(r,i),\qquad M(r_*)=D/2.
$$
Two traversals identify diameter endpoints. Distances from both endpoints then locate the centre by examining the edges. This requires $O(v)$ time and storage for $v$ stored vertices, including multifurcations. An all-zero tree has score zero and an arbitrary representative root determined by the tie rule.
Midpoint rooting uses equality of accumulated evolutionary distance as a guide to ancestry. Strong rate differences can displace the selected root even when the topology and branch lengths are known accurately. Mai, Sayyari and Mirarab give a linear-time treatment alongside their minimum-variance method [1].
## Minimum variance
For $n$ tips, define
$$
\bar d(r)=\frac{1}{n}\sum_i d(r,i),\qquad
V(r)=\frac{1}{n}\sum_i\bigl(d(r,i)-\bar d(r)\bigr)^2.
$$
Consider an edge of length $l$ directed from its stored parent $p$ to child $c$, and let $C$ be the tips on the child side. Moving the root a distance $x$ along this edge increases distances to tips outside $C$ by $x$ and decreases distances to tips inside $C$ by $x$. Write $s_i=-1$ for $i\in C$ and $s_i=1$ otherwise. Then
$$
d(r(x),i)=d(p,i)+s_i x,
$$
so
$$
V(r(x))=V(p)+B x+A x^2,
$$
where
$$
A=1-\bar s^2,\qquad
B=2\left[\frac{1}{n}\sum_i s_i d(p,i)-\bar s\,\bar d(p)\right],\qquad
\bar s=1-\frac{2|C|}{n}.
$$
For $A>0$, the minimum on this edge is
$$
x_* = \min\left(l,\max\left(0,-\frac{B}{2A}\right)\right).
$$
Subtree tip counts and distance sums provide the coefficients in constant time per edge. Propagating the quadratic to the child supplies the variance needed for the next edge. The global minimum is found in $O(v)$ time with $O(v)$ storage. This is the objective introduced by Mai, Sayyari and Mirarab [1]. Every tip has equal weight; the denominator is $n$, rather than the sample-variance denominator $n-1$.
## Minimal ancestor deviation
For a candidate root $r$, let $\alpha_r(i,j)$ be the last common ancestor of tips $i,j$ under the induced orientation. For a positive-distance pair, define
$$
\delta_{ij}(r)
=\frac{d(\alpha_r(i,j),i)-d(\alpha_r(i,j),j)}{d(i,j)}.
$$
The common path from the root to the ancestor cancels, giving the useful identity
$$
\delta_{ij}(r)=\frac{d(r,i)-d(r,j)}{d(i,j)}.
$$
The implemented MAD objective is the A1-B1 definition of Tria, Landan and Dagan [2]: the root-mean-square relative deviation over tip pairs,
$$
R(r)=\left[\frac{1}{|P|}\sum_{(i,j)\in P}\delta_{ij}(r)^2\right]^{1/2}.
$$
For strictly positive tip distances, $P$ contains all $\binom n2$ unordered pairs. On trees containing zero-length paths, Curupixa omits pairs with $d(i,j)=0$ from both the numerator and denominator. If no positive-distance pair remains, MAD returns `CRPX_ERROR_INVALID_ARGUMENT`. No small positive constant is substituted for a zero length. This convention makes the treatment of coincident tips explicit and need not match packages that regularise lengths before calculating MAD.
On an edge, pairs on the same side have constant deviations. For a pair straddling the edge, the difference of root-to-tip distances changes linearly with root position, so its squared deviation is quadratic. Summing these contributions gives
$$
S(r(x))=S(p)+B_e x+A_e x^2,\qquad S(r)=|P|R(r)^2.
$$
We minimise this quadratic on each edge, including endpoints. Its constant term includes same-side pairs as well as straddling pairs. Omitting same-side contributions would compare different objectives between candidate edges.
To accumulate coefficients efficiently, traverse from each tip $i$. For a directed edge away from $i$, let $u$ be its near endpoint, $q=d(i,u)$, and $J$ the tips beyond its far endpoint. Accumulate the subtree sums
$$
H_1=\sum_{j\in J,\ d(i,j)>0}\frac{1}{d(i,j)},\qquad
H_2=\sum_{j\in J,\ d(i,j)>0}\frac{1}{d(i,j)^2}.
$$
The contribution to the quadratic along this direction is
$$
A_i=4H_2,\qquad B_i=4(2qH_2-H_1).
$$
If this direction reverses the stored edge, substitute $y=l-x$, giving coefficients $A_i$ and $-B_i-2A_i l$ in the stored direction. Each straddling pair is encountered from both tips; dividing accumulated coefficients by two counts it once. The initial value of $S$ is calculated at the stored root and propagated across edges. Each source-tip traversal takes $O(v)$ time and reuses the same arrays, giving $O(nv)$ time and $O(v)$ scratch. For ordinary phylogenetic trees, $v=O(n)$ and the time is quadratic in tip count.
The original MAD formulation assumes contemporaneous taxa. It evaluates departures from equal ancestral distances without changing the supplied branch lengths. Applying it to serially sampled taxa would require a separate justification or an extension incorporating sampling times.
## Numerical calculations and validation
Lengths are internally divided by the maximum edge length, and accumulated quantities use `long double`. Scores for midpoint and minimum variance are converted back to their original units. A score exceeding the finite `double` range returns `CRPX_ERROR_OVERFLOW`; conversion to `double` can underflow for extremely small variance. Long double precision and range depend on the platform, and nearly tied minima should be interpreted accordingly.
`test_rooting` compares selected locations and scores with independent numerical minimisation on every edge. Its reference distances come from all-pairs shortest paths, and MAD is calculated directly from pairwise root-to-tip differences. Tests also apply the selected roots and check preservation of tip identities and pairwise distances. Fixtures cover ultrametric and non-ultrametric trees, two-tip trees, multifurcations, unary internal vertices, zero lengths and invalid inputs.
## References
[1] Mai U, Sayyari E, Mirarab S (2017). Minimum variance rooting of phylogenetic trees and implications for species tree reconstruction. *PLOS ONE* 12: e0182238. [doi:10.1371/journal.pone.0182238](https://doi.org/10.1371/journal.pone.0182238).

[2] Tria FDK, Landan G, Dagan T (2017). Phylogenetic rooting using minimal ancestor deviation. *Nature Ecology & Evolution* 1: 0193. [doi:10.1038/s41559-017-0193](https://doi.org/10.1038/s41559-017-0193). The [authors' manuscript](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad-preprint.pdf) includes the detailed algorithm and alternative definitions.
