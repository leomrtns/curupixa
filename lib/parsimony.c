/* Reimplementation of biomcmc-lib parsimony (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "parsimony.h"
#include <stdlib.h>
#include <string.h>

struct crpx_parsimony_workspace {
  crpx_tree_workspace *traversal;
  uint64_t *costs;
  size_t capacity;
};

crpx_status crpx_parsimony_workspace_new(crpx_parsimony_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_parsimony_workspace *w = calloc(1, sizeof(*w));
  if (!w) return CRPX_ERROR_ALLOCATION;
  crpx_status status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &w->traversal);
  if (status != CRPX_OK) { free(w); return status; }
  *out = w;
  return CRPX_OK;
}

void crpx_parsimony_workspace_free(crpx_parsimony_workspace *w)
{
  if (w) { crpx_tree_workspace_free(w->traversal); free(w->costs); free(w); }
}

crpx_status crpx_parsimony_score(const crpx_tree *tree, const crpx_node_id *tips, size_t tip_count,
                                unsigned int states, size_t patterns, const uint64_t *masks,
                                const size_t *weights, crpx_parsimony_workspace *w,
                                uint64_t *total, uint64_t *per_pattern)
{
  const uint64_t infinity = UINT64_MAX / 4;
  crpx_status status;
  size_t nodes = crpx_tree_node_count(tree);
  size_t *tip_index = NULL;
  uint64_t *scores = NULL;
  if (!tree || !nodes || !tips || !tip_count || !states || states > 64 || !w || !total || (!masks && patterns))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > SIZE_MAX / states || nodes * states > SIZE_MAX / sizeof(uint64_t) ||
      patterns > SIZE_MAX / tip_count || patterns > SIZE_MAX / sizeof(uint64_t) ||
      nodes > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  if (tip_count != crpx_tree_leaf_count(tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t allowed = states == 64 ? UINT64_MAX : (UINT64_C(1) << states) - 1;
  for (size_t i = 0; i < patterns * tip_count; ++i)
    if (!masks[i] || (masks[i] & ~allowed)) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_workspace_refresh(w->traversal, tree);
  if (status != CRPX_OK) { crpx_tree_workspace_reset(w->traversal); return status; }
  tip_index = malloc(nodes * sizeof(size_t));
  scores = patterns ? malloc(patterns * sizeof(uint64_t)) : NULL;
  if (!tip_index || (patterns && !scores)) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i = 0; i < nodes; ++i) tip_index[i] = SIZE_MAX;
  for (size_t i = 0; i < tip_count; ++i) {
    if (tips[i] >= nodes || crpx_tree_first_child(tree, tips[i]) != CRPX_NODE_NONE || tip_index[tips[i]] != SIZE_MAX) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
    tip_index[tips[i]] = i;
  }
  size_t required = nodes * states;
  if (required > w->capacity) {
    uint64_t *costs = realloc(w->costs, required * sizeof(uint64_t));
    if (!costs) { status = CRPX_ERROR_ALLOCATION; goto done; }
    w->costs = costs; w->capacity = required;
  }
  const crpx_node_id *order = crpx_tree_workspace_postorder(w->traversal);
  uint64_t sum = 0;
  for (size_t pattern = 0; pattern < patterns; ++pattern) {
    for (size_t i = 0; i < nodes; ++i) {
      crpx_node_id node = order[i];
      uint64_t *costs = w->costs + (size_t)node * states;
      crpx_node_id child = crpx_tree_first_child(tree, node);
      if (child == CRPX_NODE_NONE) {
        uint64_t mask = masks[pattern*tip_count+tip_index[node]];
        for (unsigned int state = 0; state < states; ++state)
          costs[state] = mask & (UINT64_C(1) << state) ? 0 : infinity;
      } else {
        for (unsigned int state = 0; state < states; ++state) costs[state] = 0;
        for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling(tree, child)) {
          uint64_t *child_costs = w->costs + (size_t)child * states, minimum = infinity;
          for (unsigned int state = 0; state < states; ++state)
            if (child_costs[state] < minimum) minimum = child_costs[state];
          for (unsigned int state = 0; state < states; ++state) {
            uint64_t addition = child_costs[state] < minimum + 1 ? child_costs[state] : minimum + 1;
            costs[state] = costs[state] >= infinity-addition ? infinity : costs[state]+addition;
          }
        }
      }
    }
    uint64_t score = infinity, *root = w->costs + (size_t)crpx_tree_root(tree)*states;
    for (unsigned int state = 0; state < states; ++state) if (root[state] < score) score = root[state];
    uint64_t weight = weights ? weights[pattern] : 1;
    if ((score && weight > UINT64_MAX/score) || sum > UINT64_MAX-score*weight) {
      status = CRPX_ERROR_OVERFLOW; goto done;
    }
    scores[pattern] = score;
    sum += score*weight;
  }
  if (per_pattern && patterns) memcpy(per_pattern, scores, patterns * sizeof(uint64_t));
  *total = sum;
done:
  crpx_tree_workspace_reset(w->traversal);
  free(tip_index); free(scores);
  return status;
}

crpx_status crpx_parsimony_score_alignment(const crpx_tree *tree, const crpx_alignment *a,
                                          crpx_parsimony_workspace *w, uint64_t *score)
{
  if (!tree || !a || !w || !score || !crpx_alignment_is_aligned(a)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa = crpx_alignment_taxa(a), patterns = crpx_alignment_patterns(a);
  if (taxa != crpx_tree_leaf_count(tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX/sizeof(crpx_node_id) || patterns > SIZE_MAX/taxa ||
      patterns*taxa > SIZE_MAX/sizeof(uint64_t)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *tips = malloc(taxa * sizeof(*tips));
  uint64_t *masks = patterns ? malloc(patterns*taxa*sizeof(*masks)) : NULL;
  bool *seen = calloc(taxa, sizeof(*seen));
  if (!tips || (patterns && !masks) || !seen) { free(tips); free(masks); free(seen); return CRPX_ERROR_ALLOCATION; }
  crpx_status status = CRPX_OK;
  for (size_t node = 0; node < crpx_tree_node_count(tree); ++node) {
    if (crpx_tree_first_child(tree, (crpx_node_id)node) != CRPX_NODE_NONE) continue;
    const char *label = crpx_tree_label(tree, (crpx_node_id)node);
    size_t taxon;
    if (!label || crpx_alignment_find_taxon(a, label, &taxon) != CRPX_OK || seen[taxon]) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
    seen[taxon] = true;
    tips[taxon] = (crpx_node_id)node;
  }
  const unsigned char *data = crpx_alignment_pattern_data(a);
  for (size_t i = 0; i < patterns*taxa; ++i) masks[i] = crpx_dna_state_mask(data[i]);
  status = crpx_parsimony_score(tree, tips, taxa, 4, patterns, masks, crpx_alignment_pattern_weights(a), w, score, NULL);
done:
  free(tips); free(masks); free(seen);
  return status;
}
