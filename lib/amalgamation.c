/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/* Szöllősi et al. (2013), Efficient Exploration of the Space of Reconciled Gene Trees, doi:10.1093/sysbio/syt054. Only
   the empirical conditional-clade component is implemented here. */
/* Build a conditional-clade distribution from sampled trees; evaluate and sample amalgamated gene trees. */
#include "amalgamation.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t left, right;
  uint64_t count;
} clade_split;
typedef struct {
  uint64_t *mask, count;
  size_t taxon, split_count;
  clade_split *splits;
} conditional_clade;
struct crpx_amalgamation {
  size_t taxa, words, clade_count, root;
  uint64_t samples;
  char **labels;
  conditional_clade *clades;
};
struct crpx_amalgamation_workspace {
  uint64_t *masks;
  size_t *clades;
  size_t node_capacity, mask_capacity;
  crpx_tree_workspace *traversal;
};

void
crpx_amalgamation_workspace_clear (crpx_amalgamation_workspace *workspace)
{
  if (!workspace) return;
  free (workspace->masks);
  free (workspace->clades);
  crpx_tree_workspace_free (workspace->traversal);
  memset (workspace, 0, sizeof (*workspace));
}

void
crpx_amalgamation_workspace_free (crpx_amalgamation_workspace *workspace)
{
  if (!workspace) return;
  crpx_amalgamation_workspace_clear (workspace);
  free (workspace);
}

crpx_status
crpx_amalgamation_workspace_new (crpx_amalgamation_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}

void
crpx_amalgamation_free (crpx_amalgamation *model)
{
  if (!model) return;
  for (size_t i = 0; i < model->clade_count; i++) {
    free (model->clades[i].mask);
    free (model->clades[i].splits);
  }
  if (model->labels)
    for (size_t i = 0; i < model->taxa; i++)
      free (model->labels[i]);
  free (model->labels);
  free (model->clades);
  free (model);
}

size_t
crpx_amalgamation_taxa (const crpx_amalgamation *model)
{
  return model ? model->taxa : 0;
}
size_t
crpx_amalgamation_clades (const crpx_amalgamation *model)
{
  return model ? model->clade_count : 0;
}
uint64_t
crpx_amalgamation_samples (const crpx_amalgamation *model)
{
  return model ? model->samples : 0;
}

static crpx_status
prepare_tree (const crpx_amalgamation *model, const crpx_tree *tree, const crpx_node_id *tips,
              crpx_amalgamation_workspace *workspace)
{
  size_t nodes = crpx_tree_node_count (tree);
  crpx_status status;
  if (!tree || !tips || crpx_tree_leaf_count (tree) != model->taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > SIZE_MAX / model->words / sizeof (uint64_t) || nodes > SIZE_MAX / sizeof (size_t))
    return CRPX_ERROR_OVERFLOW;
  if (nodes > workspace->node_capacity) {
    size_t *clades = realloc (workspace->clades, nodes * sizeof (*clades));
    if (!clades) return CRPX_ERROR_ALLOCATION;
    workspace->clades = clades;
    workspace->node_capacity = nodes;
  }
  size_t mask_words = nodes * model->words;
  if (mask_words > workspace->mask_capacity) {
    uint64_t *masks = realloc (workspace->masks, mask_words * sizeof (*masks));
    if (!masks) return CRPX_ERROR_ALLOCATION;
    workspace->masks = masks;
    workspace->mask_capacity = mask_words;
  }
  if (!workspace->traversal) {
    status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace->traversal);
    if (status != CRPX_OK) return status;
  }
  status = crpx_tree_workspace_refresh (workspace->traversal, tree);
  if (status != CRPX_OK) return status;
  memset (workspace->masks, 0, mask_words * sizeof (*workspace->masks));
  for (size_t i = 0; i < nodes; i++)
    workspace->clades[i] = SIZE_MAX;
  for (size_t taxon = 0; taxon < model->taxa; taxon++) {
    crpx_node_id node = tips[taxon];
    if (node >= nodes || crpx_tree_first_child (tree, node) != CRPX_NODE_NONE || workspace->clades[node] != SIZE_MAX)
      return CRPX_ERROR_INVALID_ARGUMENT;
    workspace->clades[node] = taxon;
    workspace->masks[(size_t) node * model->words + taxon / 64] |= UINT64_C (1) << (taxon % 64);
  }
  const crpx_node_id *order = crpx_tree_workspace_postorder (workspace->traversal);
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = order[i], left = crpx_tree_first_child (tree, node);
    if (left == CRPX_NODE_NONE) continue;
    crpx_node_id right = crpx_tree_next_sibling (tree, left);
    if (right == CRPX_NODE_NONE || crpx_tree_next_sibling (tree, right) != CRPX_NODE_NONE)
      return CRPX_ERROR_UNSUPPORTED;
    for (size_t word = 0; word < model->words; word++)
      workspace->masks[(size_t) node * model->words + word] =
        workspace->masks[(size_t) left * model->words + word] | workspace->masks[(size_t) right * model->words + word];
  }
  return CRPX_OK;
}

static size_t
find_clade (const crpx_amalgamation *model, const uint64_t *mask)
{
  for (size_t i = 0; i < model->clade_count; i++)
    if (!memcmp (model->clades[i].mask, mask, model->words * sizeof (*mask))) return i;
  return SIZE_MAX;
}

static crpx_status
record_tree (crpx_amalgamation *model, const crpx_tree *tree, crpx_amalgamation_workspace *workspace, uint64_t weight)
{
  size_t nodes = crpx_tree_node_count (tree);
  const crpx_node_id *order = crpx_tree_workspace_postorder (workspace->traversal);
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = order[i], left = crpx_tree_first_child (tree, node);
    const uint64_t *mask = workspace->masks + (size_t) node * model->words;
    size_t clade_id = find_clade (model, mask);
    if (clade_id == SIZE_MAX) {
      if (model->clade_count == SIZE_MAX / sizeof (conditional_clade)) return CRPX_ERROR_OVERFLOW;
      conditional_clade *grown = realloc (model->clades, (model->clade_count + 1) * sizeof (*grown));
      if (!grown) return CRPX_ERROR_ALLOCATION;
      model->clades = grown;
      clade_id = model->clade_count++;
      conditional_clade *clade = model->clades + clade_id;
      memset (clade, 0, sizeof (*clade));
      clade->taxon = left == CRPX_NODE_NONE ? workspace->clades[node] : SIZE_MAX;
      clade->mask = malloc (model->words * sizeof (*mask));
      if (!clade->mask) return CRPX_ERROR_ALLOCATION;
      memcpy (clade->mask, mask, model->words * sizeof (*mask));
    }
    conditional_clade *clade = model->clades + clade_id;
    clade->count += weight;
    if (left != CRPX_NODE_NONE) {
      crpx_node_id right = crpx_tree_next_sibling (tree, left);
      size_t first = workspace->clades[left], second = workspace->clades[right];
      if (first > second) {
        size_t swap = first;
        first = second;
        second = swap;
      }
      size_t split = 0;
      for (; split < clade->split_count; split++)
        if (clade->splits[split].left == first && clade->splits[split].right == second) break;
      if (split == clade->split_count) {
        if (clade->split_count == SIZE_MAX / sizeof (clade_split)) return CRPX_ERROR_OVERFLOW;
        clade_split *grown = realloc (clade->splits, (split + 1) * sizeof (*grown));
        if (!grown) return CRPX_ERROR_ALLOCATION;
        clade->splits = grown;
        clade->splits[split] = (clade_split){first, second, 0};
        clade->split_count++;
      }
      clade->splits[split].count += weight;
    }
    workspace->clades[node] = clade_id;
  }
  model->root = workspace->clades[crpx_tree_root (tree)];
  return CRPX_OK;
}

crpx_status
crpx_amalgamation_new (const crpx_tree *const *trees, const crpx_node_id *const *tip_nodes, size_t tree_count,
                       size_t taxa, const char *const *labels, const uint64_t *weights, crpx_amalgamation **out)
{ /* Count observed clades and splits once, before evaluating or sampling trees. */
  crpx_amalgamation *model;
  crpx_amalgamation_workspace workspace = {0};
  crpx_status status = CRPX_OK;
  uint64_t samples = 0;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!trees || !tip_nodes || !tree_count || !taxa || !labels) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > UINT32_MAX / 2 + (size_t) 1 || taxa > SIZE_MAX / sizeof (char *)) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < tree_count; i++) {
    uint64_t weight = weights ? weights[i] : 1;
    if (!weight) return CRPX_ERROR_INVALID_ARGUMENT;
    if (samples > UINT64_MAX - weight) return CRPX_ERROR_OVERFLOW;
    samples += weight;
  }
  for (size_t i = 0; i < taxa; i++) {
    if (!labels[i] || !*labels[i]) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t j = 0; j < i; j++)
      if (!strcmp (labels[i], labels[j])) return CRPX_ERROR_INVALID_ARGUMENT;
  }
  model = calloc (1, sizeof (*model));
  if (!model) return CRPX_ERROR_ALLOCATION;
  model->taxa = taxa;
  model->words = taxa / 64 + (taxa % 64 != 0);
  model->samples = samples;
  model->labels = calloc (taxa, sizeof (*model->labels));
  if (!model->labels) {
    crpx_amalgamation_free (model);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < taxa; i++) {
    size_t length = strlen (labels[i]) + 1;
    model->labels[i] = malloc (length);
    if (!model->labels[i]) {
      status = CRPX_ERROR_ALLOCATION;
      goto done;
    }
    memcpy (model->labels[i], labels[i], length);
  }
  for (size_t i = 0; i < tree_count; i++) {
    status = prepare_tree (model, trees[i], tip_nodes[i], &workspace);
    if (status == CRPX_OK) status = record_tree (model, trees[i], &workspace, weights ? weights[i] : 1);
    crpx_tree_workspace_reset (workspace.traversal);
    if (status != CRPX_OK) goto done;
  }
  *out = model;
  model = NULL;
done:
  crpx_amalgamation_workspace_clear (&workspace);
  crpx_amalgamation_free (model);
  return status;
}

crpx_status
crpx_amalgamation_log_probability (const crpx_amalgamation *model, const crpx_tree *tree, const crpx_node_id *tips,
                                   crpx_amalgamation_workspace *workspace, double *log_probability)
{ /* Multiply conditional split probabilities for the supplied topology in log space. */
  if (!model || !workspace || !log_probability) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare_tree (model, tree, tips, workspace);
  double result = 0;
  if (status == CRPX_OK) {
    size_t nodes = crpx_tree_node_count (tree);
    const crpx_node_id *order = crpx_tree_workspace_postorder (workspace->traversal);
    for (size_t i = 0; i < nodes; i++) {
      crpx_node_id node = order[i], left = crpx_tree_first_child (tree, node);
      size_t clade_id = find_clade (model, workspace->masks + (size_t) node * model->words);
      workspace->clades[node] = clade_id;
      if (clade_id == SIZE_MAX) {
        result = -INFINITY;
        break;
      }
      if (left == CRPX_NODE_NONE) continue;
      crpx_node_id right = crpx_tree_next_sibling (tree, left);
      size_t first = workspace->clades[left], second = workspace->clades[right];
      if (first > second) {
        size_t swap = first;
        first = second;
        second = swap;
      }
      const conditional_clade *clade = model->clades + clade_id;
      size_t split = 0;
      for (; split < clade->split_count; split++)
        if (clade->splits[split].left == first && clade->splits[split].right == second) break;
      if (split == clade->split_count) {
        result = -INFINITY;
        break;
      }
      result += log ((double) clade->splits[split].count) - log ((double) clade->count);
    }
  }
  crpx_tree_workspace_reset (workspace->traversal);
  if (status == CRPX_OK) *log_probability = result;
  return status;
}

typedef struct {
  size_t clade;
  crpx_node_id parent;
} pending_clade;

crpx_status
crpx_amalgamation_sample (const crpx_amalgamation *model, crpx_amalgamation_uniform uniform, void *state,
                          crpx_tree **out, double *log_probability)
{ /* Assemble a tree from sampled conditional splits; input sample trees need not contain this combination. */
  crpx_tree *tree = NULL;
  pending_clade *pending;
  crpx_status status;
  size_t used = 1;
  double result = 0;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!model || !uniform) return CRPX_ERROR_INVALID_ARGUMENT;
  if (model->taxa > SIZE_MAX / 2 / sizeof (*pending)) return CRPX_ERROR_OVERFLOW;
  pending = malloc (2 * model->taxa * sizeof (*pending));
  if (!pending) return CRPX_ERROR_ALLOCATION;
  pending[0] = (pending_clade){model->root, CRPX_NODE_NONE};
  status = crpx_tree_new (&tree);
  if (status == CRPX_OK) status = crpx_tree_reserve (tree, 2 * model->taxa - 1);
  while (status == CRPX_OK && used) {
    pending_clade next = pending[--used];
    const conditional_clade *clade = model->clades + next.clade;
    crpx_node_id node;
    status = crpx_tree_add_node (tree, next.parent, &node);
    if (status != CRPX_OK) break;
    if (clade->taxon != SIZE_MAX) {
      status = crpx_tree_set_label (tree, node, model->labels[clade->taxon]);
      continue;
    }
    size_t split = 0;
    if (clade->split_count > 1) {
      double draw;
      status = crpx_sample_uniform (uniform, state, &draw);
      if (status != CRPX_OK) break;
      long double target = (long double) draw * clade->count;
      uint64_t cumulative = 0;
      for (; split + 1 < clade->split_count; split++) {
        cumulative += clade->splits[split].count;
        if (target < cumulative) break;
      }
    }
    const clade_split *selected = clade->splits + split;
    result += log ((double) selected->count) - log ((double) clade->count);
    pending[used++] = (pending_clade){selected->right, node};
    pending[used++] = (pending_clade){selected->left, node};
  }
  free (pending);
  if (status != CRPX_OK) {
    crpx_tree_free (tree);
    return status;
  }
  *out = tree;
  if (log_probability) *log_probability = result;
  return CRPX_OK;
}
