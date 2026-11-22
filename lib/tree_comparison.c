/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "tree_comparison.h"
#include "tip_correspondence.h"
#include "ancestry_index.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
struct crpx_tree_comparison {
  size_t taxa, tree_tips, words, count;
  bool rooted, lengths;
  uint64_t *splits;
  double *edge_lengths, *mrca_depth, *mrca_length, *path_depth, *path_length, *pendants;
};
typedef struct {
  const uint64_t *bits;
  size_t words;
  double length;
} split_record;
static int compare_record (const void *first, const void *second);
static crpx_status prepare_pair (const crpx_tree *first, const crpx_tree *second, bool rooted, bool lengths,
                                 crpx_tree_comparison **left, crpx_tree_comparison **right);
static size_t bit_count (const uint64_t *bits, size_t words);
static bool trivial (const crpx_tree_comparison *index, size_t row);

void
crpx_tree_comparison_free (crpx_tree_comparison *index)
{
  if (!index) return;
  free (index->splits);
  free (index->edge_lengths);
  free (index->mrca_depth);
  free (index->mrca_length);
  free (index->path_depth);
  free (index->path_length);
  free (index->pendants);
  free (index);
}
crpx_status
crpx_tree_comparison_new (const crpx_tree *tree, const crpx_node_id *tips, size_t taxa, bool rooted, bool with_lengths,
                          crpx_tree_comparison **out)
{ /* Prepare coordinates once; none of the comparison functions below allocates memory. */
  crpx_tree_comparison *index = NULL;
  crpx_tree_workspace *traversal = NULL;
  crpx_ancestry_index *ancestry = NULL;
  uint64_t *descendants = NULL;
  split_record *records = NULL;
  double *heights = NULL, *lengths = NULL;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !tips || !taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t nodes = crpx_tree_node_count (tree), words = taxa / 64 + (taxa % 64 != 0);
  if (!nodes || nodes > SIZE_MAX / sizeof (split_record) || nodes > SIZE_MAX / words / sizeof (uint64_t) ||
      taxa > SIZE_MAX / taxa / sizeof (double))
    return CRPX_ERROR_OVERFLOW;
  status = crpx_ancestry_index_new (tree, &ancestry);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < taxa; i++) {
    if (!crpx_ancestry_index_is_leaf (ancestry, tips[i])) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    for (size_t j = 0; j < i; j++)
      if (tips[i] == tips[j]) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
  }
  index = calloc (1, sizeof (*index));
  if (!index) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  index->taxa = taxa;
  index->tree_tips = crpx_tree_leaf_count (tree);
  index->words = words;
  index->rooted = rooted;
  index->lengths = with_lengths;
  index->splits = calloc (nodes * words, sizeof (uint64_t));
  index->edge_lengths = calloc (nodes, sizeof (double));
  index->mrca_depth = calloc (taxa * taxa, sizeof (double));
  index->mrca_length = calloc (taxa * taxa, sizeof (double));
  index->path_depth = calloc (taxa * taxa, sizeof (double));
  index->path_length = calloc (taxa * taxa, sizeof (double));
  index->pendants = calloc (taxa, sizeof (double));
  descendants = calloc (nodes * words, sizeof (uint64_t));
  records = calloc (nodes, sizeof (*records));
  heights = calloc (nodes, sizeof (double));
  lengths = calloc (nodes, sizeof (double));
  if (!index->splits || !index->edge_lengths || !index->mrca_depth || !index->mrca_length || !index->path_depth ||
      !index->path_length || !index->pendants || !descendants || !records || !heights || !lengths) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (traversal, tree);
  if (status != CRPX_OK) goto done;
  const crpx_node_id *preorder = crpx_tree_workspace_preorder (traversal);
  const crpx_node_id *postorder = crpx_tree_workspace_postorder (traversal);
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = preorder[i], parent = crpx_tree_parent (tree, node);
    if (parent == CRPX_NODE_NONE) continue;
    if (with_lengths) {
      status = crpx_tree_branch_length (tree, node, lengths + node);
      if (status != CRPX_OK) goto done;
      if (!isfinite (lengths[node]) || lengths[node] < 0) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
    }
    heights[node] = heights[parent] + lengths[node];
    if (!isfinite (heights[node])) {
      status = CRPX_ERROR_OVERFLOW;
      goto done;
    }
  }
  for (size_t i = 0; i < taxa; i++) {
    descendants[(size_t) tips[i] * words + i / 64] |= UINT64_C (1) << (i % 64);
    index->pendants[i] = lengths[tips[i]];
    for (size_t j = 0; j <= i; j++) {
      crpx_node_id ancestor;
      status = crpx_ancestry_index_lca (ancestry, tips[i], tips[j], &ancestor);
      if (status != CRPX_OK) goto done;
      size_t coordinate = i * taxa + j;
      index->mrca_depth[coordinate] = crpx_ancestry_index_depth (ancestry, ancestor);
      index->mrca_length[coordinate] = heights[ancestor];
      index->path_depth[coordinate] = (double) crpx_ancestry_index_depth (ancestry, tips[i]) +
                                      crpx_ancestry_index_depth (ancestry, tips[j]) - 2 * index->mrca_depth[coordinate];
      if (!rooted && i != j && ancestor == crpx_tree_root (tree)) {
        crpx_node_id child = crpx_tree_first_child (tree, ancestor);
        crpx_node_id sibling = crpx_tree_next_sibling (tree, child);
        if (sibling != CRPX_NODE_NONE && crpx_tree_next_sibling (tree, sibling) == CRPX_NODE_NONE)
          index->path_depth[coordinate] -= 1;
      }
      index->path_length[coordinate] = heights[tips[i]] + heights[tips[j]] - 2 * heights[ancestor];
    }
  }
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = postorder[i], parent = crpx_tree_parent (tree, node);
    if (parent != CRPX_NODE_NONE)
      for (size_t word = 0; word < words; word++)
        descendants[(size_t) parent * words + word] |= descendants[(size_t) node * words + word];
  }
  size_t record_count = 0;
  for (size_t node = 0; node < nodes; node++) {
    if (node == crpx_tree_root (tree)) continue;
    uint64_t *bits = descendants + node * words;
    size_t ones = bit_count (bits, words);
    if (!ones || ones == taxa) continue;
    if (!rooted && (bits[0] & 1)) {
      for (size_t word = 0; word < words; word++)
        bits[word] = ~bits[word];
      if (taxa % 64) bits[words - 1] &= (UINT64_C (1) << (taxa % 64)) - 1;
    }
    records[record_count++] = (split_record){bits, words, lengths[node]};
  }
  /* Sort once, combine unary/root duplicates, then merge ordered rows during scoring. */
  qsort (records, record_count, sizeof (*records), compare_record);
  for (size_t position = 0; position < record_count; position++) {
    split_record record = records[position];
    if (!index->count || memcmp (index->splits + (index->count - 1) * words, record.bits, words * sizeof (uint64_t))) {
      memcpy (index->splits + index->count * words, record.bits, words * sizeof (uint64_t));
      index->count++;
    }
    index->edge_lengths[index->count - 1] += record.length;
  }
  *out = index;
  index = NULL;
done:
  crpx_tree_comparison_free (index);
  crpx_tree_workspace_free (traversal);
  crpx_ancestry_index_free (ancestry);
  free (descendants);
  free (records);
  free (heights);
  free (lengths);
  return status;
}
crpx_status
crpx_tree_comparison_distance (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                               crpx_tree_distance_result *result)
{
  if (!first || !second || !result || first->taxa != second->taxa || first->rooted != second->rooted)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_tree_distance_result distance = {0};
  distance.taxa = first->taxa;
  distance.first_tree_tips = first->tree_tips;
  distance.second_tree_tips = second->tree_tips;
  long double squared = 0, scale_squared = 0;
  size_t left = 0, right = 0;
  while (left < first->count || right < second->count) {
    int order = left == first->count ? 1
                : right == second->count
                  ? -1
                  : memcmp (first->splits + left * first->words, second->splits + right * second->words,
                            first->words * sizeof (uint64_t));
    double first_length = 0, second_length = 0;
    if (order <= 0) {
      bool nontrivial = !trivial (first, left);
      distance.first_splits += nontrivial;
      distance.rf += nontrivial && order < 0;
      first_length = first->edge_lengths[left++];
    }
    if (order >= 0) {
      bool nontrivial = !trivial (second, right);
      distance.second_splits += nontrivial;
      distance.rf += nontrivial && order > 0;
      second_length = second->edge_lengths[right++];
    }
    distance.length_scale += first_length + second_length;
    scale_squared += (long double) first_length * first_length + (long double) second_length * second_length;
    long double difference = (long double) first_length - second_length;
    distance.weighted_rf += (double) fabsl (difference);
    squared += difference * difference;
  }
  distance.branch_score = (double) sqrtl (squared);
  distance.branch_scale = (double) sqrtl (scale_squared);
  if (first->lengths && second->lengths &&
      (!isfinite (distance.branch_score) || !isfinite (distance.branch_scale) || !isfinite (distance.weighted_rf) ||
       !isfinite (distance.length_scale)))
    return CRPX_ERROR_OVERFLOW;
  if (!first->lengths || !second->lengths)
    distance.weighted_rf = distance.branch_score = distance.length_scale = distance.branch_scale = NAN;
  *result = distance;
  return CRPX_OK;
}
crpx_status
crpx_tree_comparison_kendall_colijn (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                                     double lambda, double *distance)
{
  if (!first || !second || !distance || !first->rooted || !second->rooted || first->taxa != second->taxa ||
      !isfinite (lambda) || lambda < 0 || lambda > 1 || (lambda && (!first->lengths || !second->lengths)))
    return CRPX_ERROR_INVALID_ARGUMENT;
  long double squared = 0;
  for (size_t i = 0; i < first->taxa; i++) {
    double terminal = lambda * (first->pendants[i] - second->pendants[i]);
    squared += (long double) terminal * terminal;
    for (size_t j = 0; j < i; j++) {
      size_t position = i * first->taxa + j;
      double difference = (1 - lambda) * (first->mrca_depth[position] - second->mrca_depth[position]) +
                          lambda * (first->mrca_length[position] - second->mrca_length[position]);
      squared += (long double) difference * difference;
    }
  }
  double value = (double) sqrtl (squared);
  if (!isfinite (value)) return CRPX_ERROR_OVERFLOW;
  *distance = value;
  return CRPX_OK;
}
crpx_status
crpx_tree_comparison_path_difference (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                                      bool weighted, double *distance)
{
  if (!first || !second || !distance || first->taxa != second->taxa || first->rooted != second->rooted ||
      (weighted && (!first->lengths || !second->lengths)))
    return CRPX_ERROR_INVALID_ARGUMENT;
  long double squared = 0;
  const double *left = weighted ? first->path_length : first->path_depth;
  const double *right = weighted ? second->path_length : second->path_depth;
  for (size_t i = 0; i < first->taxa; i++)
    for (size_t j = 0; j < i; j++) {
      long double difference = left[i * first->taxa + j] - right[i * first->taxa + j];
      squared += difference * difference;
    }
  double value = (double) sqrtl (squared);
  if (!isfinite (value)) return CRPX_ERROR_OVERFLOW;
  *distance = value;
  return CRPX_OK;
}

crpx_status
crpx_tree_weighted_distances (const crpx_tree *first, const crpx_tree *second, bool rooted,
                              crpx_tree_distance_result *result)
{
  crpx_tree_comparison *left = NULL, *right = NULL;
  crpx_status status = prepare_pair (first, second, rooted, true, &left, &right);
  if (status == CRPX_OK) status = crpx_tree_comparison_distance (left, right, result);
  crpx_tree_comparison_free (left);
  crpx_tree_comparison_free (right);
  return status;
}
crpx_status
crpx_tree_kendall_colijn_distance (const crpx_tree *first, const crpx_tree *second, double lambda, double *distance)
{
  crpx_tree_comparison *left = NULL, *right = NULL;
  crpx_status status = prepare_pair (first, second, true, lambda != 0, &left, &right);
  if (status == CRPX_OK) status = crpx_tree_comparison_kendall_colijn (left, right, lambda, distance);
  crpx_tree_comparison_free (left);
  crpx_tree_comparison_free (right);
  return status;
}
crpx_status
crpx_tree_path_difference (const crpx_tree *first, const crpx_tree *second, bool weighted, double *distance)
{
  crpx_tree_comparison *left = NULL, *right = NULL;
  crpx_status status = prepare_pair (first, second, false, weighted, &left, &right);
  if (status == CRPX_OK) status = crpx_tree_comparison_path_difference (left, right, weighted, distance);
  crpx_tree_comparison_free (left);
  crpx_tree_comparison_free (right);
  return status;
}
static crpx_status
prepare_pair (const crpx_tree *first, const crpx_tree *second, bool rooted, bool lengths, crpx_tree_comparison **left,
              crpx_tree_comparison **right)
{
  crpx_tip_correspondence *correspondence = NULL;
  crpx_status status = crpx_tip_correspondence_new (first, second, &correspondence);
  if (status == CRPX_OK)
    status = crpx_tree_comparison_new (first, correspondence->first, correspondence->taxa, rooted, lengths, left);
  if (status == CRPX_OK)
    status = crpx_tree_comparison_new (second, correspondence->second, correspondence->taxa, rooted, lengths, right);
  crpx_tip_correspondence_free (correspondence);
  return status;
}

static int
compare_record (const void *first, const void *second)
{
  const split_record *left = first, *right = second;
  return memcmp (left->bits, right->bits, left->words * sizeof (uint64_t));
}

static size_t
bit_count (const uint64_t *bits, size_t words)
{
  size_t count = 0;
  for (size_t i = 0; i < words; i++)
    for (uint64_t word = bits[i]; word; word &= word - 1)
      count++;
  return count;
}
static bool
trivial (const crpx_tree_comparison *index, size_t row)
{
  size_t count = bit_count (index->splits + row * index->words, index->words);
  return count <= 1 || (!index->rooted && count >= index->taxa - 1);
}
