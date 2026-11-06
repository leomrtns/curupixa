/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose replacements for biomcmc-lib topology randomisation.
 */
#include "topology_randomise.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t
random_bounded(crpx_tree_random_function random_word, void *state, uint64_t bound)
{
  uint64_t word, threshold = (UINT64_C(0) - bound) % bound;
  do { word = random_word(state); } while (word < threshold);
  return word % bound;
}

typedef struct { crpx_node_id first, second; } random_children;
typedef struct { crpx_node_id source, parent; } random_frame;

crpx_status
crpx_tree_random_binary(size_t leaves, const char *const *labels,
                         crpx_tree_random_function random_word, void *state, crpx_tree **result)
{
  crpx_tree *tree = NULL;
  crpx_node_id *active = NULL, node;
  random_children *children = NULL;
  random_frame *stack = NULL, frame;
  size_t count, i, first, second, pending, next;
  char label[32];
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!leaves || !random_word) return CRPX_ERROR_INVALID_ARGUMENT;
  if (leaves - 1 > (UINT32_MAX - 1) / 2 || leaves > SIZE_MAX / sizeof(*active)) return CRPX_ERROR_OVERFLOW;
  count = leaves * 2 - 1;
  if (count > SIZE_MAX / sizeof(*children) || count > SIZE_MAX / sizeof(*stack)) return CRPX_ERROR_OVERFLOW;
  if (labels) for (i = 0; i < leaves; ++i) if (!labels[i]) return CRPX_ERROR_INVALID_ARGUMENT;
  active = (crpx_node_id *) malloc(leaves * sizeof(*active));
  children = (random_children *) malloc(count * sizeof(*children));
  stack = (random_frame *) malloc(count * sizeof(*stack));
  if (!active || !children || !stack) { status = CRPX_ERROR_ALLOCATION; goto finished; }
  for (i = 0; i < leaves; ++i) active[i] = (crpx_node_id) i;
  next = leaves;
  for (i = leaves; i > 1; --i) {
    first = (size_t) random_bounded(random_word, state, i);
    node = active[first];
    active[first] = active[i - 1];
    second = (size_t) random_bounded(random_word, state, i - 1);
    children[next].first = node;
    children[next].second = active[second];
    active[second] = (crpx_node_id) next++;
  }
  status = crpx_tree_new(&tree);
  if (status != CRPX_OK) goto finished;
  status = crpx_tree_reserve(tree, count);
  if (status != CRPX_OK) goto finished;
  stack[0].source = active[0]; stack[0].parent = CRPX_NODE_NONE;
  pending = 1;
  while (pending) {
    frame = stack[--pending];
    status = crpx_tree_add_node(tree, frame.parent, &node);
    if (status != CRPX_OK) goto finished;
    if (frame.source < leaves) {
      if (!labels) snprintf(label, sizeof(label), "%u", (unsigned) frame.source);
      status = crpx_tree_set_label(tree, node, labels ? labels[frame.source] : label);
      if (status != CRPX_OK) goto finished;
    }
    else {
      stack[pending].source = children[frame.source].second; stack[pending++].parent = node;
      stack[pending].source = children[frame.source].first; stack[pending++].parent = node;
    }
  }
finished:
  free(active); free(children); free(stack);
  if (status != CRPX_OK) crpx_tree_free(tree);
  else *result = tree;
  return status;
}

static crpx_node_id
binary_sister(const crpx_tree *tree, crpx_node_id node)
{
  crpx_node_id parent = crpx_tree_parent(tree, node), first, second;
  if (parent == CRPX_NODE_NONE) return CRPX_NODE_NONE;
  first = crpx_tree_first_child(tree, parent);
  second = crpx_tree_next_sibling(tree, first);
  if (second == CRPX_NODE_NONE || crpx_tree_next_sibling(tree, second) != CRPX_NODE_NONE) return CRPX_NODE_NONE;
  return first == node ? second : first;
}

static bool
binary_node(const crpx_tree *tree, crpx_node_id node)
{
  crpx_node_id first = crpx_tree_first_child(tree, node), second;
  if (first == CRPX_NODE_NONE) return false;
  second = crpx_tree_next_sibling(tree, first);
  return second != CRPX_NODE_NONE && crpx_tree_next_sibling(tree, second) == CRPX_NODE_NONE;
}

crpx_status
crpx_tree_nni(crpx_tree *tree, crpx_node_id node, unsigned child_choice)
{
  crpx_node_id sister, child;
  if (!tree || child_choice > 1) return CRPX_ERROR_INVALID_ARGUMENT;
  if (node >= crpx_tree_node_count(tree)) return CRPX_ERROR_OUT_OF_RANGE;
  sister = binary_sister(tree, node);
  if (sister == CRPX_NODE_NONE || !binary_node(tree, node)) return CRPX_ERROR_TOPOLOGY;
  child = crpx_tree_first_child(tree, node);
  if (child_choice) child = crpx_tree_next_sibling(tree, child);
  return crpx_tree_swap_subtrees(tree, child, sister);
}

crpx_status
crpx_tree_random_nni(crpx_tree *tree, crpx_tree_random_function random_word, void *state)
{
  size_t i, eligible = 0;
  crpx_node_id selected = CRPX_NODE_NONE, node;
  if (!tree || !random_word) return CRPX_ERROR_INVALID_ARGUMENT;
  for (i = 0; i < crpx_tree_node_count(tree); ++i) {
    node = (crpx_node_id) i;
    if (!binary_node(tree, node) || binary_sister(tree, node) == CRPX_NODE_NONE) continue;
    if (!random_bounded(random_word, state, ++eligible)) selected = node;
  }
  if (selected == CRPX_NODE_NONE) return CRPX_ERROR_NOT_FOUND;
  return crpx_tree_nni(tree, selected, (unsigned) random_bounded(random_word, state, 2));
}

crpx_status
crpx_tree_random_spr(crpx_tree *tree, double graft_fraction,
                      crpx_tree_random_function random_word, void *state)
{
  crpx_tree_workspace *workspace = NULL;
  const crpx_node_id *preorder, *postorder;
  uint32_t *sizes = NULL;
  crpx_node_id selected = CRPX_NODE_NONE, node, parent, sister, target = CRPX_NODE_NONE;
  size_t count, i, eligible = 0, selected_begin = 0, selected_end, rank;
  crpx_status status;
  if (!tree || !random_word || !isfinite(graft_fraction) || graft_fraction < 0 || graft_fraction > 1)
    return CRPX_ERROR_INVALID_ARGUMENT;
  count = crpx_tree_node_count(tree);
  if (count < 4) return CRPX_ERROR_NOT_FOUND;
  sizes = (uint32_t *) malloc(count * sizeof(*sizes));
  if (!sizes) return CRPX_ERROR_ALLOCATION;
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace);
  if (status != CRPX_OK) goto finished;
  status = crpx_tree_workspace_refresh(workspace, tree);
  if (status != CRPX_OK) goto finished;
  preorder = crpx_tree_workspace_preorder(workspace);
  postorder = crpx_tree_workspace_postorder(workspace);
  for (i = 0; i < count; ++i) sizes[i] = 1;
  for (i = 0; i < count; ++i) {
    node = postorder[i];
    parent = crpx_tree_parent(tree, node);
    if (parent != CRPX_NODE_NONE) sizes[parent] += sizes[node];
  }
  for (i = 0; i < count; ++i) {
    node = preorder[i];
    if (binary_sister(tree, node) == CRPX_NODE_NONE || count - sizes[node] <= 2) continue;
    if (!random_bounded(random_word, state, ++eligible)) { selected = node; selected_begin = i; }
  }
  if (selected == CRPX_NODE_NONE) { status = CRPX_ERROR_NOT_FOUND; goto finished; }
  parent = crpx_tree_parent(tree, selected);
  sister = binary_sister(tree, selected);
  selected_end = selected_begin + sizes[selected];
  rank = (size_t) random_bounded(random_word, state, count - sizes[selected] - 2);
  for (i = 0; i < count; ++i) {
    if (i >= selected_begin && i < selected_end) continue;
    node = preorder[i];
    if (node == parent || node == sister) continue;
    if (!rank) { target = node; break; }
    --rank;
  }
  status = crpx_tree_spr(tree, selected, target, graft_fraction);
finished:
  free(sizes);
  crpx_tree_workspace_free(workspace);
  return status;
}
