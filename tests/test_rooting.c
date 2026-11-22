/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "rooting.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); exit (1); } } while (0)
#define MAX_NODES 32
static double distances[MAX_NODES][MAX_NODES];
static size_t tips[MAX_NODES], count, ntips;

static void
prepare (const crpx_tree *tree)
{ /* Independent all-pairs shortest paths, deliberately unrelated to the production traversal. */
  count = crpx_tree_node_count (tree); ntips = 0;
  CHECK (count <= MAX_NODES);
  for (size_t i = 0; i < count; i++) {
    for (size_t j = 0; j < count; j++) distances[i][j] = i == j ? 0 : INFINITY;
    if (crpx_tree_first_child (tree, i) == CRPX_NODE_NONE) tips[ntips++] = i;
  }
  /* Initialisation above must precede edge filling. */
  for (size_t i = 0; i < count; i++) {
    crpx_node_id p = crpx_tree_parent (tree, i);
    if (p != CRPX_NODE_NONE) {
      double length;
      CHECK (crpx_tree_branch_length (tree, i, &length) == CRPX_OK);
      distances[i][p] = distances[p][i] = length;
    }
  }
  for (size_t k = 0; k < count; k++) for (size_t i = 0; i < count; i++) for (size_t j = 0; j < count; j++)
    if (distances[i][k] + distances[k][j] < distances[i][j])
      distances[i][j] = distances[i][k] + distances[k][j];
}

static double
objective (const crpx_tree *tree, size_t edge, double x, unsigned method)
{
  size_t p = crpx_tree_parent (tree, edge);
  double length, d[MAX_NODES], sum = 0, result = 0, pairs = 0;
  CHECK (crpx_tree_branch_length (tree, edge, &length) == CRPX_OK);
  for (size_t i = 0; i < ntips; i++) {
    d[i] = fmin (distances[p][tips[i]] + x, distances[edge][tips[i]] + length - x);
    sum += d[i];
    if (d[i] > result) result = d[i];
  }
  if (!method) return result;
  result = 0;
  if (method == 1) {
    for (size_t i = 0; i < ntips; i++) result += (d[i] - sum / ntips) * (d[i] - sum / ntips);
    return result / ntips;
  }
  /* Shared root-to-LCA paths cancel in the difference of tip distances. */
  for (size_t i = 0; i < ntips; i++) for (size_t j = i + 1; j < ntips; j++) {
    double path = distances[tips[i]][tips[j]];
    if (path > 0) { double r = (d[i] - d[j]) / path; result += r * r; pairs++; }
  }
  return sqrt (result / pairs);
}

static void
check_tree (crpx_tree *tree, crpx_rooting_workspace *workspace)
{
  prepare (tree);
  for (unsigned method = 0; method < 3; method++) {
    crpx_root_location actual;
    uint64_t revision = crpx_tree_revision (tree);
    crpx_status status = method == 0 ? crpx_tree_root_midpoint (tree, workspace, &actual) :
                         method == 1 ? crpx_tree_root_minvar (tree, workspace, &actual) :
                                       crpx_tree_root_mad (tree, workspace, &actual);
    CHECK (status == CRPX_OK && crpx_tree_revision (tree) == revision);
    double best = INFINITY, length;
    for (size_t edge = 0; edge < count; edge++) if (edge != crpx_tree_root (tree)) {
      CHECK (crpx_tree_branch_length (tree, edge, &length) == CRPX_OK);
      double lo = 0, hi = length;
      for (unsigned iteration = 0; iteration < 100; iteration++) {
        double a = (2 * lo + hi) / 3, b = (lo + 2 * hi) / 3;
        if (objective (tree, edge, a, method) < objective (tree, edge, b, method)) hi = b;
        else lo = a;
      }
      double candidate = objective (tree, edge, (lo + hi) / 2, method);
      if (candidate < best) best = candidate;
    }
    CHECK (actual.fraction >= 0 && actual.fraction <= 1);
    CHECK (crpx_tree_branch_length (tree, actual.child, &length) == CRPX_OK);
    double measured = objective (tree, actual.child, length * actual.fraction, method);
    if (fabs (measured - best) > 2e-7 * (1 + best) || fabs (measured - actual.score) > 2e-7 * (1 + best)) {
      fprintf (stderr, "method %u: expected %.12g measured %.12g reported %.12g\n", method, best, measured, actual.score);
      exit (1);
    }
    crpx_tree *copy = NULL;
    CHECK (crpx_tree_clone (tree, &copy) == CRPX_OK);
    crpx_node_id first = crpx_tree_first_child (copy, crpx_tree_root (copy));
    bool binary = crpx_tree_next_sibling (copy, crpx_tree_next_sibling (copy, first)) == CRPX_NODE_NONE;
    crpx_node_id new_root;
    CHECK (crpx_tree_reroot_edge (copy, actual.child, actual.fraction,
             binary && ntips > 2 ? CRPX_TREE_REUSE_BINARY_ROOT : CRPX_TREE_KEEP_OLD_ROOT, &new_root) == CRPX_OK);
    CHECK (crpx_tree_validate (copy) == CRPX_OK);
    double original[MAX_NODES][MAX_NODES]; memcpy (original, distances, sizeof (original));
    size_t original_tips[MAX_NODES], original_ntips = ntips;
    memcpy (original_tips, tips, sizeof (tips));
    prepare (copy);
    CHECK (ntips == original_ntips);
    for (size_t i = 0; i < ntips; i++) for (size_t j = 0; j < ntips; j++)
      CHECK (fabs (distances[original_tips[i]][original_tips[j]] - original[original_tips[i]][original_tips[j]]) < 1e-9);
    crpx_tree_free (copy);
    prepare (tree);
  }
}

int
main (void)
{
  crpx_rooting_workspace *workspace = NULL;
  CHECK (crpx_rooting_workspace_new (&workspace) == CRPX_OK);
  const char *examples[] = {"((a:1,b:1):2,(c:2,d:2):1);", "(a:1,b:2,c:7,d:3);",
    "((a:0,b:0):2,c:3,d:1);", "(((a:1):2,b:4):1,(c:2,d:6):2);", "(a:2,b:8);",
    "((a:0,b:1):0,(c:0,d:3):0);", "(a:1,(b:2,c:3):4,(d:2,(e:1,f:6):3):1);"};
  for (size_t e = 0; e < sizeof (examples) / sizeof (*examples); e++) {
    crpx_tree *tree = NULL;
    CHECK (crpx_tree_from_newick (examples[e], &tree, NULL) == CRPX_OK);
    check_tree (tree, e % 2 ? NULL : workspace);
    if (e == 0 || e == 6) {
      for (unsigned replicate = 0; replicate < 120; replicate++) {
        for (crpx_node_id v = 0; v < crpx_tree_node_count (tree); v++) if (v != crpx_tree_root (tree))
          CHECK (crpx_tree_set_branch_length (tree, v, (double)((replicate * 71 + v * 137 + replicate * v * 17) % 101) / 13) == CRPX_OK);
        check_tree (tree, workspace);
      }
    }
    crpx_tree_free (tree);
  }
  /* Changing the storage root changes neither the objectives nor their minima. */
  crpx_tree *oriented = NULL;
  CHECK (crpx_tree_from_newick ("(a:1,(b:2,c:4):3,(d:5,e:2):1);", &oriented, NULL) == CRPX_OK);
  for (unsigned pass = 0; pass < 4; pass++) {
    check_tree (oriented, workspace);
    crpx_root_location location;
    CHECK (crpx_tree_root_minvar (oriented, workspace, &location) == CRPX_OK);
    crpx_node_id new_root;
    CHECK (crpx_tree_reroot_edge (oriented, location.child, .23,
                                 CRPX_TREE_KEEP_OLD_ROOT, &new_root) == CRPX_OK);
  }
  crpx_tree_free (oriented);
  /* Scaling all lengths preserves fractions and has the specified effect on scores. */
  CHECK (crpx_tree_from_newick ("(a:1,b:3,c:8,d:2);", &oriented, NULL) == CRPX_OK);
  crpx_root_location reference[3];
  CHECK (crpx_tree_root_midpoint (oriented, workspace, reference) == CRPX_OK);
  CHECK (crpx_tree_root_minvar (oriented, workspace, reference + 1) == CRPX_OK);
  CHECK (crpx_tree_root_mad (oriented, workspace, reference + 2) == CRPX_OK);
  for (int exponent = -100; exponent <= 100; exponent += 200) {
    crpx_tree *scaled = NULL;
    CHECK (crpx_tree_clone (oriented, &scaled) == CRPX_OK);
    double scale = pow (10, exponent);
    for (crpx_node_id v = 0; v < crpx_tree_node_count (scaled); v++) if (v != crpx_tree_root (scaled)) {
      double length;
      CHECK (crpx_tree_branch_length (scaled, v, &length) == CRPX_OK);
      CHECK (crpx_tree_set_branch_length (scaled, v, length * scale) == CRPX_OK);
    }
    crpx_root_location actual[3];
    CHECK (crpx_tree_root_midpoint (scaled, workspace, actual) == CRPX_OK);
    CHECK (crpx_tree_root_minvar (scaled, workspace, actual + 1) == CRPX_OK);
    CHECK (crpx_tree_root_mad (scaled, workspace, actual + 2) == CRPX_OK);
    for (unsigned method = 0; method < 3; method++) {
      CHECK (actual[method].child == reference[method].child);
      CHECK (fabs (actual[method].fraction - reference[method].fraction) < 1e-12);
      double factor = method == 0 ? scale : method == 1 ? scale * scale : 1;
      CHECK (fabs (actual[method].score / factor - reference[method].score) < 1e-12);
    }
    crpx_tree_free (scaled);
  }
  crpx_tree_free (oriented);
  /* Deep trees exercise iterative traversals without a quadratic distance matrix. */
  CHECK (crpx_tree_new (&oriented) == CRPX_OK);
  crpx_node_id spine;
  CHECK (crpx_tree_add_node (oriented, CRPX_NODE_NONE, &spine) == CRPX_OK);
  for (unsigned i = 0; i < 10000; i++) {
    crpx_node_id tip, next;
    CHECK (crpx_tree_add_node (oriented, spine, &tip) == CRPX_OK);
    CHECK (crpx_tree_set_branch_length (oriented, tip, 1) == CRPX_OK);
    CHECK (crpx_tree_add_node (oriented, spine, &next) == CRPX_OK);
    CHECK (crpx_tree_set_branch_length (oriented, next, 1) == CRPX_OK);
    spine = next;
  }
  crpx_root_location deep;
  CHECK (crpx_tree_root_midpoint (oriented, workspace, &deep) == CRPX_OK);
  CHECK (crpx_tree_root_minvar (oriented, workspace, &deep) == CRPX_OK);
  crpx_tree_free (oriented);
  crpx_tree *tree = NULL;
  crpx_root_location unchanged = {123, .25, 77};
  CHECK (crpx_tree_from_newick ("(a:0,b:0,c:0);", &tree, NULL) == CRPX_OK);
  CHECK (crpx_tree_root_mad (tree, workspace, &unchanged) == CRPX_ERROR_INVALID_ARGUMENT && unchanged.child == 123);
  CHECK (crpx_tree_root_midpoint (tree, workspace, &unchanged) == CRPX_OK && unchanged.score == 0);
  CHECK (crpx_tree_root_minvar (tree, workspace, &unchanged) == CRPX_OK && unchanged.score == 0);
  CHECK (crpx_tree_clear_branch_length (tree, 1) == CRPX_OK);
  CHECK (crpx_tree_root_minvar (tree, workspace, &unchanged) != CRPX_OK);
  CHECK (crpx_tree_set_branch_length (tree, 1, -1) == CRPX_OK);
  CHECK (crpx_tree_root_midpoint (tree, workspace, &unchanged) != CRPX_OK);
  crpx_tree_free (tree);
  CHECK (crpx_tree_from_newick ("((a:1,b:2):3);", &tree, NULL) == CRPX_OK);
  CHECK (crpx_tree_root_mad (tree, workspace, &unchanged) == CRPX_ERROR_UNSUPPORTED);
  crpx_tree_free (tree);
  crpx_rooting_workspace_free (workspace);
  puts ("Rooting objectives agree with independent edge searches.");
  return 0;
}
