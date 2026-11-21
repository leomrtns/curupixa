/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "dtl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* A state is an incoming species branch, or one open time segment thereof.
 * Younger states precede older states. A boundary has zero, one or two successors. */
typedef struct {
  crpx_node_id species;
  size_t first, second, interval;
} species_state;
struct crpx_dtl_species {
  size_t nodes, state_count, origin;
  species_state *states;
  size_t *subtree_first, *subtree_last, *interval_first, *interval_end;
  uint8_t *leaves;
  bool dated;
};
typedef enum { HISTORY_TIP, HISTORY_DESCEND, HISTORY_LOSS, HISTORY_SPECIATION,
               HISTORY_DUPLICATION, HISTORY_TRANSFER } history_action;
typedef struct { size_t first, second; history_action action; } history_step;
typedef struct { crpx_node_id gene; size_t species; } history_frame;
struct crpx_dtl_workspace {
  double *values;
  history_step *steps;
  crpx_node_id *first_child, *second_child;
  history_frame *stack;
  size_t cell_capacity, gene_capacity, gene_nodes;
  crpx_node_id root;
  const crpx_dtl_species *species;
  crpx_tree_workspace *traversal;
  crpx_dtl_mode mode;
  bool valid;
};

static int
compare_ages (const void *first, const void *second)
{
  double left = *(const double *) first, right = *(const double *) second;
  return (left > right) - (left < right);
}

void
crpx_dtl_species_free (crpx_dtl_species *species)
{
  if (!species) return;
  free (species->states);
  free (species->subtree_first);
  free (species->subtree_last);
  free (species->interval_first);
  free (species->interval_end);
  free (species->leaves);
  free (species);
}

crpx_status
crpx_dtl_species_new (const crpx_tree *tree, const double *ages, double origin_age, crpx_dtl_species **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!tree || !crpx_tree_node_count (tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  size_t nodes = crpx_tree_node_count (tree);
  if (nodes > (SIZE_MAX / sizeof (size_t)) - 1) return CRPX_ERROR_OVERFLOW;
  crpx_node_id root = crpx_tree_root (tree);
  if (ages && (!isfinite (origin_age) || !isfinite (ages[root]) || origin_age <= ages[root]))
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_tree_workspace *traversal = NULL;
  crpx_dtl_species *species = calloc (1, sizeof (*species));
  double *boundaries = ages ? malloc ((nodes + 1) * sizeof (*boundaries)) : NULL;
  size_t *oldest = malloc (nodes * sizeof (*oldest));
  if (!species || !oldest || (ages && !boundaries)) { status = CRPX_ERROR_ALLOCATION; goto done; }
  species->nodes = nodes;
  species->dated = ages != NULL;
  species->subtree_first = malloc (nodes * sizeof (size_t));
  species->subtree_last = malloc (nodes * sizeof (size_t));
  species->leaves = malloc (nodes * sizeof (*species->leaves));
  if (!species->subtree_first || !species->subtree_last || !species->leaves) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (traversal, tree);
  if (status != CRPX_OK) goto done;
  const crpx_node_id *preorder = crpx_tree_workspace_preorder (traversal);
  const crpx_node_id *postorder = crpx_tree_workspace_postorder (traversal);
  for (size_t position = 0; position < nodes; position++) {
    crpx_node_id node = preorder[position], parent = crpx_tree_parent (tree, node);
    species->subtree_first[node] = species->subtree_last[node] = position;
    species->leaves[node] = crpx_tree_first_child (tree, node) == CRPX_NODE_NONE;
    oldest[node] = SIZE_MAX;
    if (ages) {
      if (!isfinite (ages[node]) || ages[node] < 0 || (species->leaves[node] && ages[node] != 0) ||
          (parent != CRPX_NODE_NONE && ages[parent] <= ages[node])) {
        status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
      }
      boundaries[position] = ages[node];
    }
  }
  for (size_t position = 0; position < nodes; position++) {
    crpx_node_id node = postorder[position], parent = crpx_tree_parent (tree, node);
    if (parent != CRPX_NODE_NONE && species->subtree_last[node] > species->subtree_last[parent])
      species->subtree_last[parent] = species->subtree_last[node];
  }
  size_t intervals = 1, state_count = nodes;
  if (ages) {
    boundaries[nodes] = origin_age;
    qsort (boundaries, nodes + 1, sizeof (*boundaries), compare_ages);
    size_t distinct = 1;
    for (size_t index = 1; index <= nodes; index++)
      if (boundaries[index] != boundaries[distinct - 1]) boundaries[distinct++] = boundaries[index];
    intervals = distinct - 1;
    state_count = 0;
    for (size_t interval = 0; interval < intervals; interval++) for (size_t node = 0; node < nodes; node++) {
      crpx_node_id parent = crpx_tree_parent (tree, (crpx_node_id) node);
      double upper = parent == CRPX_NODE_NONE ? origin_age : ages[parent];
      if (ages[node] <= boundaries[interval] && upper >= boundaries[interval + 1]) {
        if (state_count == SIZE_MAX) { status = CRPX_ERROR_OVERFLOW; goto done; }
        state_count++;
      }
    }
  }
  if (state_count > SIZE_MAX / sizeof (species_state) || intervals > SIZE_MAX / sizeof (size_t)) {
    status = CRPX_ERROR_OVERFLOW; goto done;
  }
  species->states = malloc (state_count * sizeof (*species->states));
  species->interval_first = malloc (intervals * sizeof (size_t));
  species->interval_end = malloc (intervals * sizeof (size_t));
  if (!species->states || !species->interval_first || !species->interval_end) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  species->state_count = state_count;
  size_t used = 0;
  for (size_t interval = 0; interval < intervals; interval++) {
    species->interval_first[interval] = used;
    for (size_t position = 0; position < nodes; position++) {
      crpx_node_id node = ages ? (crpx_node_id) position : postorder[position];
      crpx_node_id parent = crpx_tree_parent (tree, node);
      if (ages) {
        double upper = parent == CRPX_NODE_NONE ? origin_age : ages[parent];
        if (ages[node] > boundaries[interval] || upper < boundaries[interval + 1]) continue;
      }
      species_state state = {node, SIZE_MAX, SIZE_MAX, interval};
      if (ages && ages[node] < boundaries[interval]) state.first = oldest[node];
      else {
        crpx_node_id child = crpx_tree_first_child (tree, node);
        if (child != CRPX_NODE_NONE) {
          state.first = oldest[child];
          state.second = oldest[crpx_tree_next_sibling (tree, child)];
        }
      }
      species->states[used] = state;
      oldest[node] = used++;
    }
    species->interval_end[interval] = used;
  }
  species->origin = oldest[root];
  *result = species;
  species = NULL;
done:
  free (boundaries);
  free (oldest);
  crpx_tree_workspace_free (traversal);
  crpx_dtl_species_free (species);
  return status;
}

size_t
crpx_dtl_species_states (const crpx_dtl_species *species)
{
  return species ? species->state_count : 0;
}

size_t
crpx_dtl_species_storage_size (const crpx_dtl_species *species)
{
  if (!species) return 0;
  size_t intervals = species->states[species->state_count - 1].interval + 1;
  return sizeof (*species) + species->state_count * sizeof (species_state) +
    species->nodes * (2 * sizeof (size_t) + sizeof (uint8_t)) + 2 * intervals * sizeof (size_t);
}

crpx_status
crpx_dtl_workspace_new (crpx_dtl_workspace **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = calloc (1, sizeof (**result));
  return *result ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}

void
crpx_dtl_workspace_clear (crpx_dtl_workspace *workspace)
{
  if (!workspace) return;
  free (workspace->values);
  free (workspace->steps);
  free (workspace->first_child);
  free (workspace->second_child);
  free (workspace->stack);
  crpx_tree_workspace_free (workspace->traversal);
  memset (workspace, 0, sizeof (*workspace));
}

void
crpx_dtl_workspace_free (crpx_dtl_workspace *workspace)
{
  if (!workspace) return;
  crpx_dtl_workspace_clear (workspace);
  free (workspace);
}

size_t
crpx_dtl_workspace_storage_size (const crpx_dtl_workspace *workspace)
{
  if (!workspace) return 0;
  return sizeof (*workspace) + crpx_tree_workspace_storage_size (workspace->traversal) +
    workspace->cell_capacity * (sizeof (double) + sizeof (history_step)) +
    workspace->gene_capacity * (2 * sizeof (crpx_node_id) + sizeof (history_frame));
}

static crpx_status
reserve_dtl (crpx_dtl_workspace *workspace, size_t gene_nodes, size_t species_states)
{
  if (gene_nodes > SIZE_MAX / species_states) return CRPX_ERROR_OVERFLOW;
  size_t cells = gene_nodes * species_states;
  if (cells > SIZE_MAX / sizeof (history_step) || cells > SIZE_MAX / sizeof (double) ||
      gene_nodes > SIZE_MAX / sizeof (history_frame)) return CRPX_ERROR_OVERFLOW;
  if (cells > workspace->cell_capacity) {
    double *values = malloc (cells * sizeof (*values));
    history_step *steps = malloc (cells * sizeof (*steps));
    if (!values || !steps) { free (values); free (steps); return CRPX_ERROR_ALLOCATION; }
    free (workspace->values); free (workspace->steps);
    workspace->values = values;
    workspace->steps = steps;
    workspace->cell_capacity = cells;
  }
  if (gene_nodes > workspace->gene_capacity) {
    crpx_node_id *first = malloc (gene_nodes * sizeof (*first));
    crpx_node_id *second = malloc (gene_nodes * sizeof (*second));
    history_frame *stack = malloc (gene_nodes * sizeof (*stack));
    if (!first || !second || !stack) {
      free (first); free (second); free (stack); return CRPX_ERROR_ALLOCATION;
    }
    free (workspace->first_child); free (workspace->second_child); free (workspace->stack);
    workspace->first_child = first;
    workspace->second_child = second;
    workspace->stack = stack;
    workspace->gene_capacity = gene_nodes;
  }
  return CRPX_OK;
}

/* Both modes enumerate exactly the same histories; only addition of alternative
 * histories differs. Products are sums of costs, or sums of log weights. */
static void
consider_history (crpx_dtl_mode mode, double candidate, history_step step, double *value, history_step *chosen)
{
  if (mode == CRPX_DTL_MINIMUM_COST) {
    if (candidate < *value) { *value = candidate; *chosen = step; }
  } else if (candidate != -INFINITY) {
    if (*value == -INFINITY) *value = candidate;
    else {
      double larger = *value > candidate ? *value : candidate;
      double smaller = *value > candidate ? candidate : *value;
      *value = larger + log1p (exp (smaller - larger));
    }
  }
}

static bool
transfer_allowed (const crpx_dtl_species *species, size_t donor, size_t recipient)
{
  if (donor == recipient) return false;
  if (species->dated) return true;
  crpx_node_id first = species->states[donor].species, second = species->states[recipient].species;
  bool second_below_first = species->subtree_first[first] <= species->subtree_first[second] &&
    species->subtree_first[second] <= species->subtree_last[first];
  bool first_below_second = species->subtree_first[second] <= species->subtree_first[first] &&
    species->subtree_first[first] <= species->subtree_last[second];
  return !second_below_first && !first_below_second;
}

crpx_status
crpx_dtl_evaluate (const crpx_tree *gene, const crpx_dtl_species *species, const crpx_node_id *tip_map,
                   const crpx_dtl_costs *costs, crpx_dtl_mode mode, crpx_dtl_workspace *workspace, double *value)
{
  if (workspace) workspace->valid = false;
  if (!gene || !species || !tip_map || !costs || !workspace || !value ||
      (mode != CRPX_DTL_MINIMUM_COST && mode != CRPX_DTL_LOG_PARTITION) ||
      !isfinite (costs->duplication) || costs->duplication < 0 ||
      !isfinite (costs->transfer) || costs->transfer < 0 || !isfinite (costs->loss) || costs->loss < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (gene, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  size_t gene_nodes = crpx_tree_node_count (gene), states = species->state_count;
  if (!gene_nodes) return CRPX_ERROR_INVALID_ARGUMENT;
  status = reserve_dtl (workspace, gene_nodes, states);
  if (status != CRPX_OK) return status;
  if (!workspace->traversal)
    status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace->traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (workspace->traversal, gene);
  if (status != CRPX_OK) { crpx_tree_workspace_reset (workspace->traversal); return status; }
  const crpx_node_id *postorder = crpx_tree_workspace_postorder (workspace->traversal);
  double duplication_weight = mode == CRPX_DTL_MINIMUM_COST ? costs->duplication : -costs->duplication;
  double transfer_weight = mode == CRPX_DTL_MINIMUM_COST ? costs->transfer : -costs->transfer;
  double loss_weight = mode == CRPX_DTL_MINIMUM_COST ? costs->loss : -costs->loss;
  for (size_t position = 0; position < gene_nodes; position++) {
    crpx_node_id node = postorder[position], first = crpx_tree_first_child (gene, node);
    crpx_node_id second = first == CRPX_NODE_NONE ? CRPX_NODE_NONE : crpx_tree_next_sibling (gene, first);
    workspace->first_child[node] = first;
    workspace->second_child[node] = second;
    if (first == CRPX_NODE_NONE && (tip_map[node] >= species->nodes || !species->leaves[tip_map[node]])) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
    double *row = workspace->values + node * states;
    const double *left = first == CRPX_NODE_NONE ? NULL : workspace->values + first * states;
    const double *right = second == CRPX_NODE_NONE ? NULL : workspace->values + second * states;
    for (size_t state_id = 0; state_id < states; state_id++) {
      species_state state = species->states[state_id];
      double best = mode == CRPX_DTL_MINIMUM_COST ? INFINITY : -INFINITY;
      history_step chosen = {SIZE_MAX, SIZE_MAX, HISTORY_TIP};
      if (state.first == SIZE_MAX) {
        if (first == CRPX_NODE_NONE && tip_map[node] == state.species) best = 0;
      } else if (state.second == SIZE_MAX) {
        consider_history (mode, row[state.first], (history_step) {state.first, SIZE_MAX, HISTORY_DESCEND},
                          &best, &chosen);
      } else {
        consider_history (mode, loss_weight + row[state.first],
                          (history_step) {state.first, SIZE_MAX, HISTORY_LOSS}, &best, &chosen);
        consider_history (mode, loss_weight + row[state.second],
                          (history_step) {state.second, SIZE_MAX, HISTORY_LOSS}, &best, &chosen);
        if (first != CRPX_NODE_NONE) {
          consider_history (mode, left[state.first] + right[state.second],
                            (history_step) {state.first, state.second, HISTORY_SPECIATION}, &best, &chosen);
          consider_history (mode, left[state.second] + right[state.first],
                            (history_step) {state.second, state.first, HISTORY_SPECIATION}, &best, &chosen);
        }
      }
      if (first != CRPX_NODE_NONE) {
        consider_history (mode, duplication_weight + left[state_id] + right[state_id],
                          (history_step) {state_id, state_id, HISTORY_DUPLICATION}, &best, &chosen);
        for (size_t recipient = species->interval_first[state.interval];
             recipient < species->interval_end[state.interval]; recipient++) {
          if (!transfer_allowed (species, state_id, recipient)) continue;
          consider_history (mode, transfer_weight + left[state_id] + right[recipient],
                            (history_step) {state_id, recipient, HISTORY_TRANSFER}, &best, &chosen);
          consider_history (mode, transfer_weight + left[recipient] + right[state_id],
                            (history_step) {recipient, state_id, HISTORY_TRANSFER}, &best, &chosen);
        }
      }
      row[state_id] = best;
      workspace->steps[node * states + state_id] = chosen;
    }
  }
  {
    crpx_node_id root = crpx_tree_root (gene);
    double result = workspace->values[root * states + species->origin];
    if (!isfinite (result)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    workspace->root = root;
    workspace->gene_nodes = gene_nodes;
    workspace->species = species;
    workspace->mode = mode;
    workspace->valid = true;
    *value = result;
  }
done:
  crpx_tree_workspace_reset (workspace->traversal);
  return status;
}

crpx_status
crpx_dtl_traceback (crpx_dtl_workspace *workspace, crpx_dtl_mapping *mapping,
                    size_t capacity, crpx_dtl_counts *counts)
{
  if (!workspace || !workspace->valid || workspace->mode != CRPX_DTL_MINIMUM_COST || !mapping || !counts)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity < workspace->gene_nodes) return CRPX_ERROR_OUT_OF_RANGE;
  const crpx_dtl_species *species = workspace->species;
  crpx_dtl_counts total = {0, 0, 0};
  size_t used = 1;
  workspace->stack[0] = (history_frame) {workspace->root, species->origin};
  while (used) {
    history_frame frame = workspace->stack[--used];
    history_step step = workspace->steps[frame.gene * species->state_count + frame.species];
    while (step.action == HISTORY_DESCEND || step.action == HISTORY_LOSS) {
      if (step.action == HISTORY_LOSS) total.losses++;
      frame.species = step.first;
      step = workspace->steps[frame.gene * species->state_count + frame.species];
    }
    crpx_dtl_mapping entry = {CRPX_DTL_TIP, species->states[frame.species].species, CRPX_NODE_NONE};
    if (step.action != HISTORY_TIP) {
      if (step.action == HISTORY_SPECIATION) entry.event = CRPX_DTL_SPECIATION;
      else if (step.action == HISTORY_DUPLICATION) { entry.event = CRPX_DTL_DUPLICATION; total.duplications++; }
      else {
        entry.event = CRPX_DTL_TRANSFER;
        total.transfers++;
        size_t recipient = step.first == frame.species ? step.second : step.first;
        entry.recipient = species->states[recipient].species;
      }
      workspace->stack[used++] = (history_frame) {workspace->first_child[frame.gene], step.first};
      workspace->stack[used++] = (history_frame) {workspace->second_child[frame.gene], step.second};
    }
    mapping[frame.gene] = entry;
  }
  *counts = total;
  return CRPX_OK;
}
