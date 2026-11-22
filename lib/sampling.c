/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "sampling.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static crpx_status bounded (size_t bound, crpx_uniform_function uniform, void *state, size_t *value);

crpx_status
crpx_sample_uniform (crpx_uniform_function uniform, void *state, double *value)
{
  if (!uniform || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  double draw = uniform (state);
  if (!isfinite (draw) || draw < 0 || draw >= 1) return CRPX_ERROR_INVALID_ARGUMENT;
  *value = draw;
  return CRPX_OK;
}

crpx_status
crpx_sample_log_weights (const double *weights, size_t count, crpx_uniform_function uniform, void *state,
                         size_t *selected)
{ /* Rescale before exponentiation so very small conditional probabilities remain usable. */
  if (!weights || !count || !selected) return CRPX_ERROR_INVALID_ARGUMENT;
  double maximum = -INFINITY, draw;
  long double total = 0, cumulative = 0;
  for (size_t i = 0; i < count; i++) {
    if (isnan (weights[i]) || weights[i] == INFINITY) return CRPX_ERROR_INVALID_ARGUMENT;
    if (weights[i] > maximum) maximum = weights[i];
  }
  if (!isfinite (maximum)) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_sample_uniform (uniform, state, &draw);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < count; i++)
    total += expl ((long double) weights[i] - maximum);
  size_t last = 0;
  for (size_t i = 0; i < count; i++)
    if (isfinite (weights[i])) {
      last = i;
      cumulative += expl ((long double) weights[i] - maximum);
      if ((long double) draw * total < cumulative) {
        *selected = i;
        return CRPX_OK;
      }
    }
  *selected = last;
  return CRPX_OK;
}

crpx_status
crpx_sample_permutation (size_t *indices, size_t count, crpx_uniform_function uniform, void *state)
{
  if ((count && !indices) || !uniform) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < count; i++)
    indices[i] = i;
  for (size_t i = count; i > 1; i--) {
    size_t index;
    crpx_status status = bounded (i, uniform, state, &index);
    if (status != CRPX_OK) return status;
    size_t temporary = indices[index];
    indices[index] = indices[i - 1];
    indices[i - 1] = temporary;
  }
  return CRPX_OK;
}

typedef struct {
  size_t source;
  crpx_node_id parent;
} history_frame;
crpx_status
crpx_tree_from_binary_history (size_t leaves, const size_t *children, size_t root, const double *ages,
                               const char *const *labels, crpx_tree **out)
{
  crpx_tree *tree = NULL;
  history_frame *stack = NULL;
  unsigned char *seen = NULL;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!leaves || (leaves > 1 && !children)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (leaves > UINT32_MAX / 2 || leaves > SIZE_MAX / 2 / sizeof (*stack)) return CRPX_ERROR_OVERFLOW;
  size_t nodes = 2 * leaves - 1;
  if (root >= nodes) return CRPX_ERROR_INVALID_ARGUMENT;
  stack = malloc (nodes * sizeof (*stack));
  seen = calloc (nodes, 1);
  if (!stack || !seen) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  status = crpx_tree_new (&tree);
  size_t pending = 1, visited = 0;
  stack[0] = (history_frame){root, CRPX_NODE_NONE};
  while (status == CRPX_OK && pending) {
    history_frame frame = stack[--pending];
    if (frame.source >= nodes || seen[frame.source]) {
      status = CRPX_ERROR_TOPOLOGY;
      break;
    }
    seen[frame.source] = 1;
    visited++;
    crpx_node_id node;
    status = crpx_tree_add_node (tree, frame.parent, &node);
    if (status != CRPX_OK) break;
    if (ages) {
      if (!isfinite (ages[frame.source]) || ages[frame.source] < 0) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        break;
      }
    }
    if (frame.source < leaves) {
      char number[32];
      snprintf (number, sizeof (number), "%zu", frame.source);
      status = crpx_tree_set_label (tree, node, labels ? labels[frame.source] : number);
    } else {
      size_t first = children[2 * (frame.source - leaves)], second = children[2 * (frame.source - leaves) + 1];
      if (pending > nodes - 2 || first >= nodes || second >= nodes) {
        status = CRPX_ERROR_TOPOLOGY;
        break;
      }
      stack[pending++] = (history_frame){second, node};
      stack[pending++] = (history_frame){first, node};
    }
  }
  if (status == CRPX_OK && visited != nodes) status = CRPX_ERROR_TOPOLOGY;
  /* Traversal assigned IDs in preorder; replay it to attach edge lengths. */
  if (status == CRPX_OK && ages) {
    size_t pending = 1, next = 0;
    stack[0] = (history_frame){root, CRPX_NODE_NONE};
    while (pending && status == CRPX_OK) {
      history_frame frame = stack[--pending];
      crpx_node_id node = (crpx_node_id) next++;
      if (frame.parent != CRPX_NODE_NONE) {
        double length = ages[frame.parent] - ages[frame.source];
        if (!isfinite (length) || length < 0) {
          status = CRPX_ERROR_INVALID_ARGUMENT;
          break;
        }
        status = crpx_tree_set_branch_length (tree, node, length);
      }
      if (frame.source >= leaves) {
        stack[pending++] = (history_frame){children[2 * (frame.source - leaves) + 1], (crpx_node_id) frame.source};
        stack[pending++] = (history_frame){children[2 * (frame.source - leaves)], (crpx_node_id) frame.source};
      }
    }
  }
done:
  free (stack);
  free (seen);
  if (status == CRPX_OK)
    *out = tree;
  else
    crpx_tree_free (tree);
  return status;
}

crpx_status
crpx_tree_resolve_random (const crpx_tree *source, crpx_uniform_function uniform, void *state, crpx_tree **out)
{ /* Resolve with local edge insertions in arrays, then construct the tree once. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!source || !uniform) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate (source);
  if (status != CRPX_OK) return status;
  size_t original = crpx_tree_node_count (source);
  if (!original) return CRPX_ERROR_INVALID_ARGUMENT;
  if (original > UINT32_MAX / 2 || original > SIZE_MAX / 2 / sizeof (double) ||
      original > SIZE_MAX / 2 / sizeof (const char *))
    return CRPX_ERROR_OVERFLOW;
  size_t capacity = 2 * original, used_nodes = original;
  crpx_node_id *parents = malloc (capacity * sizeof (*parents));
  crpx_node_id *children = malloc (capacity * sizeof (*children));
  crpx_node_id *local = malloc (capacity * sizeof (*local));
  const char **labels = calloc (capacity, sizeof (*labels));
  double *lengths = calloc (capacity, sizeof (*lengths));
  if (!parents || !children || !local || !labels || !lengths) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t node = 0; node < original; node++) {
    parents[node] = crpx_tree_parent (source, (crpx_node_id) node);
    labels[node] = crpx_tree_label (source, (crpx_node_id) node);
    status = crpx_tree_branch_length (source, (crpx_node_id) node, lengths + node);
    if (status == CRPX_ERROR_NOT_FOUND) {
      lengths[node] = NAN;
      status = CRPX_OK;
    }
    if (status != CRPX_OK) goto done;
  }
  for (size_t index = 0; index < original; index++) {
    crpx_node_id node = (crpx_node_id) index;
    size_t count = 0;
    for (crpx_node_id child = crpx_tree_first_child (source, node); child != CRPX_NODE_NONE;
         child = crpx_tree_next_sibling (source, child))
      children[count++] = child;
    if (count <= 2) continue;
    crpx_node_id root_first = children[0], root_second = children[1];
    size_t used = 3;
    local[0] = node;
    local[1] = root_first;
    local[2] = root_second;
    for (size_t tip = 2; tip < count; tip++) {
      size_t selected;
      status = bounded (used, uniform, state, &selected);
      if (status != CRPX_OK) goto done;
      crpx_node_id edge = local[selected], inserted = (crpx_node_id) used_nodes++;
      if (edge == node) {
        parents[inserted] = node;
        parents[root_first] = parents[root_second] = inserted;
        root_first = inserted;
        root_second = children[tip];
      } else {
        parents[inserted] = parents[edge];
        parents[edge] = inserted;
        parents[children[tip]] = inserted;
        if (edge == root_first) root_first = inserted;
        if (edge == root_second) root_second = inserted;
      }
      local[used++] = inserted;
      local[used++] = children[tip];
    }
  }
  status = crpx_tree_from_parents (parents, used_nodes, labels, lengths, out);
done:
  free (parents);
  free (children);
  free (local);
  free (labels);
  free (lengths);
  return status;
}

crpx_status
crpx_tree_from_parents (const crpx_node_id *parents, size_t nodes, const char *const *labels, const double *lengths,
                        crpx_tree **out)
{ /* Counting offsets avoid repeated scans when constructing wide polytomies. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!nodes) return crpx_tree_new (out);
  if (!parents) return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > UINT32_MAX || nodes > SIZE_MAX / sizeof (history_frame) || nodes >= SIZE_MAX / sizeof (size_t))
    return CRPX_ERROR_OVERFLOW;
  size_t *offsets = calloc (nodes + 1, sizeof (*offsets)), *positions = malloc (nodes * sizeof (*positions));
  crpx_node_id *children = malloc (nodes * sizeof (*children));
  history_frame *stack = malloc (nodes * sizeof (*stack));
  crpx_tree *tree = NULL;
  crpx_tree_builder *builder = NULL;
  crpx_status status = CRPX_OK;
  if (!offsets || !positions || !children || !stack) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  crpx_node_id root = CRPX_NODE_NONE;
  for (size_t node = 0; node < nodes; node++) {
    if (lengths && !isnan (lengths[node]) && !isfinite (lengths[node])) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if (parents[node] == CRPX_NODE_NONE) {
      if (root != CRPX_NODE_NONE) {
        status = CRPX_ERROR_TOPOLOGY;
        goto done;
      }
      root = (crpx_node_id) node;
    } else {
      if (parents[node] >= nodes || parents[node] == node) {
        status = CRPX_ERROR_TOPOLOGY;
        goto done;
      }
      offsets[parents[node] + 1]++;
    }
  }
  if (root == CRPX_NODE_NONE) {
    status = CRPX_ERROR_TOPOLOGY;
    goto done;
  }
  for (size_t node = 1; node <= nodes; node++)
    offsets[node] += offsets[node - 1];
  memcpy (positions, offsets, nodes * sizeof (*positions));
  for (size_t node = 0; node < nodes; node++)
    if (parents[node] != CRPX_NODE_NONE) children[positions[parents[node]]++] = (crpx_node_id) node;
  status = crpx_tree_new (&tree);
  if (status == CRPX_OK) status = crpx_tree_builder_new (&builder);
  size_t pending = 1, visited = 0;
  stack[0] = (history_frame){root, CRPX_NODE_NONE};
  while (pending && status == CRPX_OK) {
    history_frame frame = stack[--pending];
    crpx_node_id node;
    status = crpx_tree_builder_append (builder, tree, frame.parent, &node);
    if (status != CRPX_OK) break;
    visited++;
    if (labels && labels[frame.source]) status = crpx_tree_set_label (tree, node, labels[frame.source]);
    if (status == CRPX_OK && lengths && !isnan (lengths[frame.source]))
      status = crpx_tree_set_branch_length (tree, node, lengths[frame.source]);
    for (size_t position = offsets[frame.source + 1]; position > offsets[frame.source];)
      stack[pending++] = (history_frame){children[--position], node};
  }
  if (status == CRPX_OK && visited != nodes) status = CRPX_ERROR_TOPOLOGY;
  if (status == CRPX_OK) {
    *out = tree;
    tree = NULL;
  }
done:
  crpx_tree_builder_free (builder);
  crpx_tree_free (tree);
  free (offsets);
  free (positions);
  free (children);
  free (stack);
  return status;
}

static crpx_status
bounded (size_t bound, crpx_uniform_function uniform, void *state, size_t *value)
{
  double draw;
  crpx_status status = crpx_sample_uniform (uniform, state, &draw);
  if (status == CRPX_OK) *value = (size_t) ((long double) draw * bound);
  return status;
}

crpx_status
crpx_tree_insert_tip (const crpx_tree *source, crpx_node_id child, double fraction, double pendant_length,
                      double stem_length, const char *label, crpx_tree **out)
{ /* Build the augmented parent array once, including the above-root alternative. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t nodes = crpx_tree_node_count (source);
  if (!source || child >= nodes || !label || !*label || !isfinite (fraction) || fraction <= 0 || fraction >= 1 ||
      !isfinite (pendant_length) || pendant_length < 0 || !isfinite (stem_length) || stem_length < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > UINT32_MAX - 2 || nodes > SIZE_MAX / sizeof (const char *) - 2 || nodes > SIZE_MAX / sizeof (double) - 2)
    return CRPX_ERROR_OVERFLOW;
  crpx_status status =
    crpx_tree_validate_policy (source, CRPX_TREE_REQUIRE_COMPLETE_LENGTHS | CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS);
  if (status != CRPX_OK) return status;
  crpx_node_id *parents = malloc ((nodes + 2) * sizeof (*parents));
  const char **labels = calloc (nodes + 2, sizeof (*labels));
  double *lengths = calloc (nodes + 2, sizeof (*lengths));
  if (!parents || !labels || !lengths) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t node = 0; node < nodes; node++) {
    parents[node] = crpx_tree_parent (source, (crpx_node_id) node);
    labels[node] = crpx_tree_label (source, (crpx_node_id) node);
    if (parents[node] != CRPX_NODE_NONE) (void) crpx_tree_branch_length (source, (crpx_node_id) node, lengths + node);
  }
  parents[nodes] = parents[child];
  parents[child] = (crpx_node_id) nodes;
  parents[nodes + 1] = (crpx_node_id) nodes;
  labels[nodes + 1] = label;
  lengths[nodes + 1] = pendant_length;
  if (parents[nodes] == CRPX_NODE_NONE)
    lengths[child] = stem_length;
  else {
    lengths[nodes] = lengths[child] * fraction;
    lengths[child] *= 1 - fraction;
  }
  status = crpx_tree_from_parents (parents, nodes + 2, labels, lengths, out);
done:
  free (parents);
  free (labels);
  free (lengths);
  return status;
}
