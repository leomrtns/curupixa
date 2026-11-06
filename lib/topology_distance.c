/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose replacements for biomcmc-lib tree distance operations.
 */
#include "topology_distance.h"
#include "hashtable.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

crpx_status
crpx_tree_is_ancestor(const crpx_tree *tree, crpx_node_id ancestor, crpx_node_id node, bool *result)
{
  if (!tree || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  if (ancestor >= crpx_tree_node_count(tree) || node >= crpx_tree_node_count(tree)) return CRPX_ERROR_OUT_OF_RANGE;
  while (node != CRPX_NODE_NONE && node != ancestor) node = crpx_tree_parent(tree, node);
  *result = node == ancestor;
  return CRPX_OK;
}

crpx_status
crpx_tree_lca(const crpx_tree *tree, crpx_node_id first, crpx_node_id second, crpx_node_id *result)
{
  size_t first_depth = 0, second_depth = 0;
  crpx_node_id node;
  if (!tree || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  if (first >= crpx_tree_node_count(tree) || second >= crpx_tree_node_count(tree)) return CRPX_ERROR_OUT_OF_RANGE;
  for (node = first; node != CRPX_NODE_NONE; node = crpx_tree_parent(tree, node)) ++first_depth;
  for (node = second; node != CRPX_NODE_NONE; node = crpx_tree_parent(tree, node)) ++second_depth;
  while (first_depth > second_depth) { first = crpx_tree_parent(tree, first); --first_depth; }
  while (second_depth > first_depth) { second = crpx_tree_parent(tree, second); --second_depth; }
  while (first != second) { first = crpx_tree_parent(tree, first); second = crpx_tree_parent(tree, second); }
  *result = first;
  return CRPX_OK;
}

static crpx_status
edge_length(const crpx_tree *tree, crpx_node_id node, double missing_length, double *result)
{
  crpx_status status = crpx_tree_branch_length(tree, node, result);
  if (status == CRPX_ERROR_NOT_FOUND && isfinite(missing_length)) { *result = missing_length; return CRPX_OK; }
  return status;
}

crpx_status
crpx_tree_patristic_distance(const crpx_tree *tree, crpx_node_id first, crpx_node_id second,
                             double missing_length, double *result)
{
  crpx_node_id lca, node;
  double length, sum = 0;
  crpx_status status;
  if (!tree || !result || isinf(missing_length)) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_lca(tree, first, second, &lca);
  if (status != CRPX_OK) return status;
  for (node = first; node != lca; node = crpx_tree_parent(tree, node)) {
    status = edge_length(tree, node, missing_length, &length);
    if (status != CRPX_OK) return status;
    sum += length;
    if (!isfinite(sum)) return CRPX_ERROR_OVERFLOW;
  }
  for (node = second; node != lca; node = crpx_tree_parent(tree, node)) {
    status = edge_length(tree, node, missing_length, &length);
    if (status != CRPX_OK) return status;
    sum += length;
    if (!isfinite(sum)) return CRPX_ERROR_OVERFLOW;
  }
  *result = sum;
  return CRPX_OK;
}

typedef struct {
  crpx_node_id node, from;
  double distance;
} distance_frame;

crpx_status
crpx_tree_patristic_matrix(const crpx_tree *tree, double missing_length,
                           crpx_distance_matrix_t *result, crpx_node_id **leaf_order)
{
  crpx_distance_matrix_t matrix = NULL;
  crpx_node_id *leaves = NULL, *indices = NULL, node, child, parent;
  distance_frame *stack = NULL, current;
  double *lengths = NULL, length;
  size_t count, leaf_count, i, j, pending;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (leaf_order) *leaf_order = NULL;
  if (!tree || isinf(missing_length) || (!isnan(missing_length) && missing_length < 0))
    return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  count = crpx_tree_node_count(tree);
  if (!count) return CRPX_ERROR_TOPOLOGY;
  leaf_count = crpx_tree_leaf_count(tree);
  if (count > SIZE_MAX / sizeof(*stack) || count > SIZE_MAX / sizeof(*lengths) ||
      count > SIZE_MAX / sizeof(*indices) || leaf_count > SIZE_MAX / sizeof(*leaves)) return CRPX_ERROR_OVERFLOW;
  stack = (distance_frame *) malloc(count * sizeof(*stack));
  lengths = (double *) malloc(count * sizeof(*lengths));
  indices = (crpx_node_id *) malloc(count * sizeof(*indices));
  leaves = (crpx_node_id *) malloc(leaf_count * sizeof(*leaves));
  if (!stack || !lengths || !indices || !leaves) { status = CRPX_ERROR_ALLOCATION; goto finished; }
  j = 0;
  for (i = 0; i < count; ++i) {
    node = (crpx_node_id) i;
    indices[i] = CRPX_NODE_NONE;
    if (crpx_tree_first_child(tree, node) == CRPX_NODE_NONE) {
      leaves[j] = node;
      indices[i] = (crpx_node_id) j++;
    }
    if (node == crpx_tree_root(tree)) lengths[i] = 0;
    else {
      status = edge_length(tree, node, missing_length, &lengths[i]);
      if (status != CRPX_OK) goto finished;
      if (lengths[i] < 0) { status = CRPX_ERROR_INVALID_ARGUMENT; goto finished; }
    }
  }
  status = crpx_distance_matrix_new(leaf_count, &matrix);
  if (status != CRPX_OK) goto finished;
  for (i = 0; i < leaf_count; ++i) {
    stack[0].node = leaves[i]; stack[0].from = CRPX_NODE_NONE; stack[0].distance = 0;
    pending = 1;
    while (pending) {
      current = stack[--pending];
      node = current.node;
      if (indices[node] != CRPX_NODE_NONE && indices[node] > i) {
        status = crpx_distance_matrix_set(matrix, i, indices[node], current.distance);
        if (status != CRPX_OK) goto finished;
      }
      parent = crpx_tree_parent(tree, node);
      if (parent != CRPX_NODE_NONE && parent != current.from) {
        length = current.distance + lengths[node];
        if (!isfinite(length)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
        stack[pending].node = parent; stack[pending].from = node; stack[pending++].distance = length;
      }
      for (child = crpx_tree_first_child(tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling(tree, child)) {
        if (child == current.from) continue;
        length = current.distance + lengths[child];
        if (!isfinite(length)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
        stack[pending].node = child; stack[pending].from = node; stack[pending++].distance = length;
      }
    }
  }
finished:
  free(stack); free(lengths); free(indices);
  if (status == CRPX_OK) {
    *result = matrix;
    if (leaf_order) *leaf_order = leaves;
    else free(leaves);
  }
  else { crpx_distance_matrix_free(matrix); free(leaves); }
  return status;
}

void crpx_tree_node_ids_free(crpx_node_id *nodes) { free(nodes); }

static size_t
word_popcount(uint64_t word)
{
  size_t count = 0;
  while (word) { word &= word - 1; ++count; }
  return count;
}

static crpx_status
collect_splits(const crpx_tree *tree, crpx_hashtable_t labels, size_t leaf_count, bool rooted,
               bool second, crpx_hashtable_t splits, size_t *distance)
{
  uint64_t *bits = NULL, *node_bits, *parent_bits, label_id, state;
  uint8_t *seen = NULL;
  crpx_tree_workspace *workspace = NULL;
  const crpx_node_id *postorder;
  crpx_node_id node, parent;
  size_t count = crpx_tree_node_count(tree), words = leaf_count / 64 + (leaf_count % 64 != 0);
  size_t i, j, ones;
  uint64_t mask = leaf_count % 64 ? (UINT64_C(1) << (leaf_count % 64)) - 1 : UINT64_MAX;
  crpx_status status;
  const char *label;
  bool replaced;
  if (count > SIZE_MAX / words || count * words > SIZE_MAX / sizeof(*bits)) return CRPX_ERROR_OVERFLOW;
  bits = (uint64_t *) calloc(count * words, sizeof(*bits));
  seen = (uint8_t *) calloc(leaf_count, sizeof(*seen));
  if (!bits || !seen) { status = CRPX_ERROR_ALLOCATION; goto finished; }
  for (i = 0; i < count; ++i) {
    if (crpx_tree_first_child(tree, (crpx_node_id) i) != CRPX_NODE_NONE) continue;
    label = crpx_tree_label(tree, (crpx_node_id) i);
    if (!label) { status = CRPX_ERROR_INVALID_ARGUMENT; goto finished; }
    status = crpx_hashtable_get_string(labels, label, &label_id);
    if (status == CRPX_ERROR_NOT_FOUND) status = CRPX_ERROR_INVALID_ARGUMENT;
    if (status != CRPX_OK) goto finished;
    if (seen[label_id]) { status = CRPX_ERROR_INVALID_ARGUMENT; goto finished; }
    seen[label_id] = 1;
    bits[i * words + label_id / 64] |= UINT64_C(1) << (label_id % 64);
  }
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace);
  if (status != CRPX_OK) goto finished;
  status = crpx_tree_workspace_refresh(workspace, tree);
  if (status != CRPX_OK) goto finished;
  postorder = crpx_tree_workspace_postorder(workspace);
  for (i = 0; i < count; ++i) {
    node = postorder[i];
    parent = crpx_tree_parent(tree, node);
    if (parent == CRPX_NODE_NONE) continue;
    node_bits = bits + node * words;
    parent_bits = bits + parent * words;
    for (j = 0; j < words; ++j) parent_bits[j] |= node_bits[j];
  }
  for (i = 0; i < count; ++i) {
    node_bits = bits + i * words;
    ones = 0;
    for (j = 0; j < words; ++j) ones += word_popcount(node_bits[j]);
    if (ones <= 1 || ones == leaf_count || (!rooted && ones >= leaf_count - 1)) continue;
    if (!rooted && (node_bits[0] & 1)) {
      for (j = 0; j < words; ++j) node_bits[j] = ~node_bits[j];
      node_bits[words - 1] &= mask;
    }
    if (!second) {
      status = crpx_hashtable_put(splits, node_bits, words * sizeof(*bits), 1, &replaced);
      if (status != CRPX_OK) goto finished;
      if (!replaced) ++*distance;
    }
    else {
      status = crpx_hashtable_get(splits, node_bits, words * sizeof(*bits), &state);
      if (status == CRPX_ERROR_NOT_FOUND) {
        if (*distance == SIZE_MAX) { status = CRPX_ERROR_OVERFLOW; goto finished; }
        status = crpx_hashtable_put(splits, node_bits, words * sizeof(*bits), 2, NULL);
        if (status != CRPX_OK) goto finished;
        ++*distance;
      }
      else if (status != CRPX_OK) goto finished;
      else if (state == 1) {
        status = crpx_hashtable_put(splits, node_bits, words * sizeof(*bits), 0, NULL);
        if (status != CRPX_OK) goto finished;
        --*distance;
      }
    }
  }
  status = CRPX_OK;
finished:
  crpx_tree_workspace_free(workspace);
  free(bits); free(seen);
  return status;
}

crpx_status
crpx_tree_rf_distance(const crpx_tree *first, const crpx_tree *second, bool rooted, size_t *result)
{
  crpx_hashtable_t labels = NULL, splits = NULL;
  size_t leaves, i, label_id = 0, distance = 0;
  const char *label;
  bool replaced;
  crpx_status status;
  if (!first || !second || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_validate(first);
  if (status != CRPX_OK) return status;
  status = crpx_tree_validate(second);
  if (status != CRPX_OK) return status;
  leaves = crpx_tree_leaf_count(first);
  if (!leaves || leaves != crpx_tree_leaf_count(second)) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_hashtable_new(leaves, 0, &labels);
  if (status != CRPX_OK) goto finished;
  for (i = 0; i < crpx_tree_node_count(first); ++i) {
    if (crpx_tree_first_child(first, (crpx_node_id) i) != CRPX_NODE_NONE) continue;
    label = crpx_tree_label(first, (crpx_node_id) i);
    if (!label) { status = CRPX_ERROR_INVALID_ARGUMENT; goto finished; }
    status = crpx_hashtable_put_string(labels, label, label_id++, &replaced);
    if (status != CRPX_OK) goto finished;
    if (replaced) { status = CRPX_ERROR_INVALID_ARGUMENT; goto finished; }
  }
  status = crpx_hashtable_new(crpx_tree_node_count(first), 0, &splits);
  if (status != CRPX_OK) goto finished;
  status = collect_splits(first, labels, leaves, rooted, false, splits, &distance);
  if (status != CRPX_OK) goto finished;
  status = collect_splits(second, labels, leaves, rooted, true, splits, &distance);
  if (status == CRPX_OK) *result = distance;
finished:
  crpx_hashtable_free(labels);
  crpx_hashtable_free(splits);
  return status;
}
