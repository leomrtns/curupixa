/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Search tree topologies by NNI using an objective supplied by the application. */
#include "tree_search.h"
#include <math.h>

static crpx_status
evaluate (const crpx_tree *tree, crpx_tree_objective objective, void *state,
          double *score, size_t *evaluations)
{
  crpx_status status;
  if (*evaluations == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  status = objective (tree, state, score);
  (*evaluations)++;
  if (status == CRPX_OK && !isfinite (*score)) return CRPX_ERROR_INVALID_ARGUMENT;
  return status;
}

static crpx_status
search_nni (const crpx_tree *start, crpx_tree_objective objective, void *state,
            size_t max_moves, bool suppress_root, crpx_tree_search_result *result)
{ /* Evaluate neighbouring topologies with the supplied objective and accept only strict improvements. */
  crpx_tree_search_result search = {0};
  crpx_status status;
  if (!start || !objective || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_validate_policy (start, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  status = crpx_tree_clone (start, &search.tree);
  if (status != CRPX_OK) return status;
  status = evaluate (search.tree, objective, state, &search.score, &search.evaluations);
  if (status != CRPX_OK) goto failed;
  while (!max_moves || search.moves < max_moves) {
    double best_score = search.score;
    crpx_node_id best_child = CRPX_NODE_NONE, best_sibling = CRPX_NODE_NONE;
    for (size_t index = 0; index < crpx_tree_node_count (search.tree); index++) {
      crpx_node_id node = (crpx_node_id) index, child = crpx_tree_first_child (search.tree, node);
      crpx_node_id parent = crpx_tree_parent (search.tree, node), sibling;
      if (parent == CRPX_NODE_NONE || child == CRPX_NODE_NONE) continue;
      sibling = crpx_tree_first_child (search.tree, parent);
      if (sibling == node) sibling = crpx_tree_next_sibling (search.tree, sibling);
      crpx_node_id children[] = {child, crpx_tree_next_sibling (search.tree, child)};
      for (size_t side = 0; side < 2; side++) {
        double score;
        status = crpx_tree_swap_subtrees (search.tree, children[side], sibling);
        if (status != CRPX_OK) goto failed;
        status = evaluate (search.tree, objective, state, &score, &search.evaluations);
        if (status != CRPX_OK) goto failed;
        status = crpx_tree_swap_subtrees (search.tree, children[side], sibling);
        if (status != CRPX_OK) goto failed;
        if (score < best_score) {
          best_score = score;
          best_child = children[side];
          best_sibling = sibling;
        }
      }
    }
    /* Suppress the artificial binary root conceptually. Its two incident edges form one unrooted edge, whose NNI
   alternatives exchange grandchildren on opposite sides. Ordinary rooted moves only relocate this root. */
    if (suppress_root) {
      crpx_node_id left = crpx_tree_first_child (search.tree, crpx_tree_root (search.tree));
      crpx_node_id right = crpx_tree_next_sibling (search.tree, left);
      if (left != CRPX_NODE_NONE && right != CRPX_NODE_NONE) {
        crpx_node_id a = crpx_tree_first_child (search.tree, left);
        crpx_node_id b = crpx_tree_first_child (search.tree, right);
        if (a != CRPX_NODE_NONE && b != CRPX_NODE_NONE) {
          crpx_node_id alternatives[] = {b, crpx_tree_next_sibling (search.tree, b)};
          for (size_t side = 0; side < 2; side++) {
            double score;
            status = crpx_tree_swap_subtrees (search.tree, a, alternatives[side]);
            if (status != CRPX_OK) goto failed;
            status = evaluate (search.tree, objective, state, &score, &search.evaluations);
            if (status != CRPX_OK) goto failed;
            status = crpx_tree_swap_subtrees (search.tree, a, alternatives[side]);
            if (status != CRPX_OK) goto failed;
            if (score < best_score) {
              best_score = score; best_child = a; best_sibling = alternatives[side];
            }
          }
        }
      }
    }
    if (best_child == CRPX_NODE_NONE) { search.local_optimum = true; break; }
    status = crpx_tree_swap_subtrees (search.tree, best_child, best_sibling);
    if (status != CRPX_OK) goto failed;
    search.score = best_score;
    search.moves++;
  }
  /* Leave callback-owned numerical caches consistent with the returned tree, rather than the last rejected neighbour. */
  status = evaluate (search.tree, objective, state, &search.score, &search.evaluations);
  if (status != CRPX_OK) goto failed;
  *result = search;
  return CRPX_OK;
failed:
  crpx_tree_free (search.tree);
  return status;
}

crpx_status
crpx_tree_search_nni (const crpx_tree *start, crpx_tree_objective objective,
                      void *state, size_t max_moves, crpx_tree_search_result *result)
{
  return search_nni (start, objective, state, max_moves, false, result);
}

crpx_status
crpx_tree_search_unrooted_nni (const crpx_tree *start, crpx_tree_objective objective,
                               void *state, size_t max_moves, crpx_tree_search_result *result)
{
  return search_nni (start, objective, state, max_moves, true, result);
}
