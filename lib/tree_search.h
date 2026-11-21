/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_TREE_SEARCH_H
#define CURUPIXA_TREE_SEARCH_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef crpx_status (*crpx_tree_objective) (const crpx_tree *tree, void *state, double *score);
typedef struct {
  crpx_tree *tree;
  double score;
  size_t evaluations, moves;
  bool local_optimum;
} crpx_tree_search_result;
/* Minimise a caller-supplied finite score by steepest descent over rooted binary
 * NNI neighbours. Each move swaps one child of a non-root internal node with its
 * sibling. Branch lengths remain attached to node IDs, as in swap_subtrees().
 * Starting tree is cloned once, never modified. Stable tip IDs preserve explicit
 * mappings; callback state/workspaces remain caller-owned and may be reused.
 * Callback must be deterministic for a tree and its fixed parameters, must not
 * mutate/retain the trial tree or cache solely by its address. The final call
 * evaluates the returned topology, so callback-owned caches describe that tree.
 * max_moves=0 continues until no strictly improving neighbour exists; otherwise
 * the move limit can stop before local optimality is established. This is a
 * generic local optimiser, not the constrained dynamic programme of ASTRAL.
 * On success result owns tree (release with crpx_tree_free); on failure result
 * is unchanged. No file handling, random draws or internal worker threads. */
crpx_status crpx_tree_search_nni (const crpx_tree *start, crpx_tree_objective objective,
                                 void *state, size_t max_moves, crpx_tree_search_result *result);
/* Same search with NNI alternatives across the edge represented by the two
 * root-child edges. Intended for root-invariant objectives such as unrooted RF
 * and quartet support. Input storage is still binary with a degree-two root;
 * its root location has no biological meaning. No lengths are combined. */
crpx_status crpx_tree_search_unrooted_nni (const crpx_tree *start, crpx_tree_objective objective,
                                          void *state, size_t max_moves, crpx_tree_search_result *result);
#ifdef __cplusplus
}
#endif
#endif
