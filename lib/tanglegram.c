/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/* Minimise crossings between explicit tip links by reordering both trees, including unequal tree sizes. */
#include "tanglegram.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t nodes, leaves;
  crpx_node_id root;
  size_t *offset;
  crpx_node_id *children, *initial, *stack, *order, *best_order;
  size_t *rank;
} ordered_tree;

struct crpx_tanglegram_workspace {
  ordered_tree first, second;
  size_t *link_offsets, *fenwick;
  crpx_node_id *targets;
  size_t count, evaluations, limit;
  uint64_t best;
};

static void
free_ordered_tree (ordered_tree *tree)
{
  free (tree->offset); free (tree->children); free (tree->initial);
  free (tree->stack); free (tree->order); free (tree->best_order); free (tree->rank);
}

void
crpx_tanglegram_workspace_free (crpx_tanglegram_workspace *workspace)
{
  if (!workspace) return;
  free_ordered_tree (&workspace->first);
  free_ordered_tree (&workspace->second);
  free (workspace->link_offsets);
  free (workspace->fenwick);
  free (workspace->targets);
  free (workspace);
}

static crpx_status
snapshot_order (const crpx_tree *source, ordered_tree *tree)
{
  crpx_status status = crpx_tree_validate (source);
  size_t used = 0;
  if (status != CRPX_OK) return status;
  tree->nodes = crpx_tree_node_count (source);
  tree->leaves = crpx_tree_leaf_count (source);
  tree->root = crpx_tree_root (source);
  if (!tree->nodes) return CRPX_ERROR_TOPOLOGY;
  if (tree->nodes > SIZE_MAX / sizeof (size_t) - 1 ||
      tree->nodes > SIZE_MAX / sizeof (crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  tree->offset = malloc ((tree->nodes + 1) * sizeof (*tree->offset));
  tree->children = malloc (tree->nodes * sizeof (*tree->children));
  tree->initial = malloc (tree->nodes * sizeof (*tree->initial));
  tree->stack = malloc (tree->nodes * sizeof (*tree->stack));
  tree->order = malloc (tree->leaves * sizeof (*tree->order));
  tree->best_order = malloc (tree->leaves * sizeof (*tree->best_order));
  tree->rank = malloc (tree->nodes * sizeof (*tree->rank));
  if (!tree->offset || !tree->children || !tree->initial || !tree->stack ||
      !tree->order || !tree->best_order || !tree->rank) return CRPX_ERROR_ALLOCATION;
  for (size_t node = 0; node < tree->nodes; node++) {
    tree->offset[node] = used;
    for (crpx_node_id child = crpx_tree_first_child (source, (crpx_node_id) node);
         child != CRPX_NODE_NONE; child = crpx_tree_next_sibling (source, child)) tree->children[used++] = child;
  }
  tree->offset[tree->nodes] = used;
  memcpy (tree->initial, tree->children, used * sizeof (*tree->initial));
  return CRPX_OK;
}

crpx_status
crpx_tanglegram_workspace_new (const crpx_tree *first, const crpx_tree *second,
                              const crpx_tanglegram_link *links, size_t count,
                              crpx_tanglegram_workspace **result)
{ /* Snapshot both topologies and group links by endpoint for repeated crossing counts. */
  crpx_tanglegram_workspace *workspace;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!first || !second || (count && !links)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / sizeof (*links)) return CRPX_ERROR_OVERFLOW;
  /* This conservative bound also prevents the crossing counter overflowing. */
  if (count > UINT32_MAX) return CRPX_ERROR_OVERFLOW;
  workspace = calloc (1, sizeof (*workspace));
  if (!workspace) return CRPX_ERROR_ALLOCATION;
  status = snapshot_order (first, &workspace->first);
  if (status == CRPX_OK) status = snapshot_order (second, &workspace->second);
  if (status != CRPX_OK) goto failed;
  for (size_t i = 0; i < count; i++) {
    if (links[i].first >= workspace->first.nodes || links[i].second >= workspace->second.nodes ||
        crpx_tree_first_child (first, links[i].first) != CRPX_NODE_NONE ||
        crpx_tree_first_child (second, links[i].second) != CRPX_NODE_NONE) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto failed;
    }
  }
  workspace->link_offsets = calloc (workspace->first.nodes + 1, sizeof (size_t));
  workspace->fenwick = calloc (workspace->second.leaves + 1, sizeof (size_t));
  workspace->targets = count ? malloc (count * sizeof (*workspace->targets)) : NULL;
  if (!workspace->link_offsets || !workspace->fenwick || (count && !workspace->targets)) {
    status = CRPX_ERROR_ALLOCATION;
    goto failed;
  }
  for (size_t i = 0; i < count; i++) workspace->link_offsets[links[i].first + 1]++;
  for (size_t node = 1; node <= workspace->first.nodes; node++)
    workspace->link_offsets[node] += workspace->link_offsets[node - 1];
  /* The traversal stack is free scratch during construction. Fill backwards so the offsets themselves remain available
   without another allocation. */
  for (size_t node = 0; node < workspace->first.nodes; node++)
    workspace->first.stack[node] = (crpx_node_id) workspace->link_offsets[node + 1];
  for (size_t i = 0; i < count; i++)
    workspace->targets[--workspace->first.stack[links[i].first]] = links[i].second;
  workspace->count = count;
  *result = workspace;
  return CRPX_OK;
failed:
  crpx_tanglegram_workspace_free (workspace);
  return status;
}

static void
leaf_order (ordered_tree *tree)
{
  size_t pending = 1, leaves = 0;
  tree->stack[0] = tree->root;
  while (pending) {
    crpx_node_id node = tree->stack[--pending];
    size_t start = tree->offset[node], end = tree->offset[node + 1];
    if (start == end) { tree->rank[node] = leaves; tree->order[leaves++] = node; }
    else while (end > start) tree->stack[pending++] = tree->children[--end];
  }
}

static uint64_t
evaluate_order (crpx_tanglegram_workspace *workspace)
{
  uint64_t crossings = 0;
  leaf_order (&workspace->first);
  leaf_order (&workspace->second);
  size_t previous = 0;
  memset (workspace->fenwick, 0, (workspace->second.leaves + 1) * sizeof (size_t));
  for (size_t leaf = 0; leaf < workspace->first.leaves; leaf++) {
    crpx_node_id node = workspace->first.order[leaf];
    size_t start = workspace->link_offsets[node], end = workspace->link_offsets[node + 1];
    for (size_t link = start; link < end; link++) {
      size_t rank = workspace->second.rank[workspace->targets[link]] + 1, prefix = 0;
      for (size_t position = rank; position; position &= position - 1) prefix += workspace->fenwick[position];
      crossings += previous - prefix;
    }
    /* Query the whole group before inserting it: links sharing this endpoint cannot cross. Inclusive prefix sums similarly
   exclude equal right tips. */
    for (size_t link = start; link < end; link++) {
      size_t rank = workspace->second.rank[workspace->targets[link]] + 1;
      for (size_t position = rank; position <= workspace->second.leaves; position += position & (~position + 1))
        workspace->fenwick[position]++;
    }
    previous += end - start;
  }
  workspace->evaluations++;
  if (crossings < workspace->best) {
    workspace->best = crossings;
    memcpy (workspace->first.best_order, workspace->first.order,
            workspace->first.leaves * sizeof (crpx_node_id));
    memcpy (workspace->second.best_order, workspace->second.order,
            workspace->second.leaves * sizeof (crpx_node_id));
  }
  return crossings;
}

static bool
budget_available (const crpx_tanglegram_workspace *workspace)
{
  return workspace->evaluations < SIZE_MAX && (!workspace->limit || workspace->evaluations < workspace->limit);
}

static void
swap_children (crpx_node_id *first, crpx_node_id *second)
{
  crpx_node_id temporary = *first;
  *first = *second;
  *second = temporary;
}

static void
improve_orders (crpx_tanglegram_workspace *workspace, uint64_t current)
{
  ordered_tree *trees[] = {&workspace->first, &workspace->second};
  bool improved = true;
  while (improved && workspace->best && budget_available (workspace)) {
    improved = false;
    for (size_t side = 0; side < 2; side++) {
      ordered_tree *tree = trees[side];
      for (size_t node = 0; node < tree->nodes; node++)
        for (size_t i = tree->offset[node]; i < tree->offset[node + 1]; i++)
          for (size_t j = i + 1; j < tree->offset[node + 1]; j++) {
            if (!workspace->best || !budget_available (workspace)) return;
            swap_children (tree->children + i, tree->children + j);
            uint64_t candidate = evaluate_order (workspace);
            if (candidate < current) { current = candidate; improved = true; }
            else swap_children (tree->children + i, tree->children + j);
          }
    }
  }
}

static int
compare_nodes (const void *first, const void *second)
{
  crpx_node_id a = *(const crpx_node_id *) first, b = *(const crpx_node_id *) second;
  return (a > b) - (a < b);
}

/* Advance one digit of a mixed-radix enumeration of child permutations. Exhausted digits return to ascending order
   before carrying to the next. */
static bool
next_order (ordered_tree *tree)
{
  for (size_t node = 0; node < tree->nodes; node++) {
    size_t start = tree->offset[node], end = tree->offset[node + 1], pivot = end;
    if (end - start < 2) continue;
    pivot = end - 1;
    while (pivot > start && tree->children[pivot - 1] >= tree->children[pivot]) pivot--;
    if (pivot > start) {
      size_t successor = end - 1;
      while (tree->children[successor] <= tree->children[pivot - 1]) successor--;
      swap_children (tree->children + pivot - 1, tree->children + successor);
    }
    size_t left = pivot, right = end;
    while (left < right && left < --right) swap_children (tree->children + left++, tree->children + right);
    if (pivot > start) return true;
  }
  return false;
}

crpx_status
crpx_tanglegram_optimise (crpx_tanglegram_workspace *workspace, size_t max_evaluations,
                         crpx_tanglegram_result *result)
{ /* Improve both orders, then enumerate permutations until completion or the evaluation limit. */
  ordered_tree *trees[2];
  bool complete = false;
  if (!workspace || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  trees[0] = &workspace->first;
  trees[1] = &workspace->second;
  workspace->evaluations = 0;
  workspace->limit = max_evaluations;
  workspace->best = UINT64_MAX;
  for (size_t side = 0; side < 2; side++)
    memcpy (trees[side]->children, trees[side]->initial,
            trees[side]->offset[trees[side]->nodes] * sizeof (crpx_node_id));
  improve_orders (workspace, evaluate_order (workspace));
  if (workspace->best && budget_available (workspace)) {
    for (size_t side = 0; side < 2; side++)
      for (size_t node = 0; node < trees[side]->nodes; node++)
        qsort (trees[side]->children + trees[side]->offset[node],
               trees[side]->offset[node + 1] - trees[side]->offset[node], sizeof (crpx_node_id), compare_nodes);
    do {
      evaluate_order (workspace);
      complete = !next_order (trees[0]) && !next_order (trees[1]);
    } while (!complete && workspace->best && budget_available (workspace));
  }
  *result = (crpx_tanglegram_result) {workspace->best, workspace->evaluations,
    complete || workspace->best == 0, trees[0]->best_order, trees[1]->best_order,
    trees[0]->leaves, trees[1]->leaves};
  return CRPX_OK;
}
