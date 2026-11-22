/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "phylogenetic_regression.h"
#include "ancestry_index.h"
#include <stdlib.h>
#include <math.h>

typedef struct {
  const crpx_tree *tree;
  crpx_ancestry_index *ancestry;
  double *heights;
  size_t nodes;
  crpx_trait_covariance_options options;
} trait_context;
static crpx_status prepare_context (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                                    const crpx_trait_covariance_options *options, trait_context *context);
static void clear_context (trait_context *context);
static double covariance_entry (const trait_context *context, crpx_node_id first, crpx_node_id second);
static crpx_status observed_covariance (const trait_context *context, const crpx_node_id *tips, size_t count,
                                        crpx_covariance **out);

void
crpx_trait_covariance_options_init (crpx_trait_covariance_options *options)
{
  if (options) *options = (crpx_trait_covariance_options){CRPX_TRAIT_BROWNIAN, 1, 1, NULL};
}

crpx_status
crpx_phylogenetic_covariance (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                              const crpx_trait_covariance_options *options, crpx_covariance **out)
{ /* Separate tree/model preparation from the general regression calculation. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  trait_context context = {0};
  crpx_status status = prepare_context (tree, tips, count, options, &context);
  if (status == CRPX_OK) status = observed_covariance (&context, tips, count, out);
  clear_context (&context);
  return status;
}

crpx_status
crpx_phylogenetic_prediction (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                              const crpx_trait_covariance_options *options, crpx_gaussian_prediction **out)
{ /* The same covariance entries condition latent node values on noisy observations. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  trait_context context = {0};
  crpx_covariance *observed = NULL;
  double *cross = NULL, *variances = NULL;
  crpx_status status = prepare_context (tree, tips, count, options, &context);
  if (status != CRPX_OK) goto done;
  if (context.nodes > SIZE_MAX / count / sizeof (double)) {
    status = CRPX_ERROR_OVERFLOW;
    goto done;
  }
  status = observed_covariance (&context, tips, count, &observed);
  if (status != CRPX_OK) goto done;
  cross = malloc (context.nodes * count * sizeof (double));
  variances = malloc (context.nodes * sizeof (double));
  if (!cross || !variances) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t node = 0; node < context.nodes; node++) {
    variances[node] = covariance_entry (&context, (crpx_node_id) node, (crpx_node_id) node);
    for (size_t tip = 0; tip < count; tip++)
      cross[node * count + tip] = covariance_entry (&context, (crpx_node_id) node, tips[tip]);
  }
  status = crpx_gaussian_prediction_new (observed, cross, variances, context.nodes, out);
done:
  crpx_covariance_free (observed);
  clear_context (&context);
  free (cross);
  free (variances);
  return status;
}

static void
clear_context (trait_context *context)
{
  crpx_ancestry_index_free (context->ancestry);
  free (context->heights);
}

static crpx_status
prepare_context (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                 const crpx_trait_covariance_options *options, trait_context *context)
{
  crpx_trait_covariance_options_init (&context->options);
  if (options) context->options = *options;
  options = &context->options;
  if (!tree || !tips || !count || !isfinite (options->lambda) || options->lambda < 0 || options->lambda > 1 ||
      options->kind < CRPX_TRAIT_BROWNIAN || options->kind > CRPX_TRAIT_OU_STATIONARY ||
      (options->kind == CRPX_TRAIT_OU_STATIONARY && (!isfinite (options->alpha) || options->alpha <= 0)))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / count / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  context->tree = tree;
  context->nodes = crpx_tree_node_count (tree);
  crpx_status status = crpx_ancestry_index_new (tree, &context->ancestry);
  if (status != CRPX_OK) return status;
  context->heights = calloc (context->nodes, sizeof (double));
  if (!context->heights) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < count; i++) {
    if (!crpx_ancestry_index_is_leaf (context->ancestry, tips[i])) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t j = 0; j < i; j++)
      if (tips[i] == tips[j]) return CRPX_ERROR_INVALID_ARGUMENT;
    if (options->tip_error_variances &&
        (!isfinite (options->tip_error_variances[i]) || options->tip_error_variances[i] < 0))
      return CRPX_ERROR_INVALID_ARGUMENT;
  }
  crpx_tree_workspace *traversal = NULL;
  status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (traversal, tree);
  if (status == CRPX_OK) {
    const crpx_node_id *order = crpx_tree_workspace_preorder (traversal);
    for (size_t i = 0; i < context->nodes; i++) {
      crpx_node_id node = order[i], parent = crpx_tree_parent (tree, node);
      if (parent == CRPX_NODE_NONE) continue;
      double length;
      status = crpx_tree_branch_length (tree, node, &length);
      if (status != CRPX_OK) break;
      if (!isfinite (length) || length < 0) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        break;
      }
      context->heights[node] = context->heights[parent] + length;
      if (!isfinite (context->heights[node])) {
        status = CRPX_ERROR_OVERFLOW;
        break;
      }
    }
  }
  crpx_tree_workspace_free (traversal);
  return status;
}

static double
covariance_entry (const trait_context *context, crpx_node_id first, crpx_node_id second)
{
  crpx_node_id ancestor;
  /* Indices were validated during preparation; every query addresses the snapshot. */
  (void) crpx_ancestry_index_lca (context->ancestry, first, second, &ancestor);
  double value = context->heights[ancestor];
  if (context->options.kind == CRPX_TRAIT_OU_STATIONARY) {
    double distance =
      (context->heights[first] - context->heights[ancestor]) + (context->heights[second] - context->heights[ancestor]);
    value = exp (-context->options.alpha * distance) / (2 * context->options.alpha);
  }
  return first == second ? value : value * context->options.lambda;
}

static crpx_status
observed_covariance (const trait_context *context, const crpx_node_id *tips, size_t count, crpx_covariance **out)
{
  double *matrix = malloc (count * count * sizeof (double));
  if (!matrix) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < count; i++)
    for (size_t j = 0; j <= i; j++) {
      double value = covariance_entry (context, tips[i], tips[j]);
      if (i == j && context->options.tip_error_variances) value += context->options.tip_error_variances[i];
      matrix[i * count + j] = matrix[j * count + i] = value;
    }
  crpx_status status = crpx_covariance_new (matrix, count, out);
  free (matrix);
  return status;
}
