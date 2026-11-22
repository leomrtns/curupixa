/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "molecular_clock.h"
#include "prob_distribution.h"
#include <math.h>
static crpx_status prepare (const crpx_clock_model *model, const crpx_tree *tree, const double *ages,
                            crpx_tree_workspace *workspace, size_t count);
static double variance (const crpx_clock_model *model, const crpx_tree *tree, const double *ages, crpx_node_id node);
crpx_status
crpx_clock_sample (const crpx_clock_model *model, const crpx_tree *tree, const double *ages,
                   crpx_uniform_function uniform, void *state, crpx_tree_workspace *workspace, double *rates,
                   size_t count)
{ /* Traverse ancestors first so autocorrelated draws use the sampled parent rate. */
  if (!rates || rates == ages) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare (model, tree, ages, workspace, count);
  if (status != CRPX_OK) return status;
  const crpx_node_id *order = crpx_tree_workspace_preorder (workspace);
  for (size_t i = 0; i < count; i++) {
    crpx_node_id node = order[i], parent = crpx_tree_parent (tree, node);
    double mean =
      parent == CRPX_NODE_NONE || model->kind != CRPX_CLOCK_AUTOCORRELATED ? model->mean_rate : rates[parent];
    double v = parent == CRPX_NODE_NONE ? 0 : variance (model, tree, ages, node);
    if (parent == CRPX_NODE_NONE || model->kind == CRPX_CLOCK_STRICT ||
        (model->kind == CRPX_CLOCK_AUTOCORRELATED && !v))
      rates[node] = mean;
    else {
      double u;
      status = crpx_sample_uniform (uniform, state, &u);
      if (status != CRPX_OK) break;
      if (!u) u = 0x1p-54;
      rates[node] = model->kind == CRPX_CLOCK_EXPONENTIAL ? -mean * log1p (-u)
                                                          : crpx_qlnorm (u, log (mean) - v / 2, sqrt (v), false);
      if (!isfinite (rates[node]) || rates[node] <= 0) {
        status = CRPX_ERROR_OVERFLOW;
        break;
      }
    }
  }
  crpx_tree_workspace_reset (workspace);
  return status;
}
crpx_status
crpx_clock_log_density (const crpx_clock_model *model, const crpx_tree *tree, const double *ages, const double *rates,
                        crpx_tree_workspace *workspace, size_t count, double *value)
{ /* The same conditional factors define sampling and scoring. */
  if (!rates || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare (model, tree, ages, workspace, count);
  if (status != CRPX_OK) return status;
  double total = 0;
  for (size_t i = 0; i < count; i++) {
    crpx_node_id node = (crpx_node_id) i, parent = crpx_tree_parent (tree, node);
    if (!isfinite (rates[node]) || rates[node] <= 0) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      break;
    }
    double mean =
      parent == CRPX_NODE_NONE || model->kind != CRPX_CLOCK_AUTOCORRELATED ? model->mean_rate : rates[parent];
    double v = parent == CRPX_NODE_NONE ? 0 : variance (model, tree, ages, node);
    if (parent == CRPX_NODE_NONE || model->kind == CRPX_CLOCK_STRICT ||
        (model->kind == CRPX_CLOCK_AUTOCORRELATED && !v)) {
      if (rates[node] != mean) total = -INFINITY;
    } else if (model->kind == CRPX_CLOCK_EXPONENTIAL)
      total += -log (mean) - rates[node] / mean;
    else
      total += crpx_dlnorm (rates[node], log (mean) - v / 2, sqrt (v), true);
  }
  crpx_tree_workspace_reset (workspace);
  if (status != CRPX_OK) return status;
  if (isnan (total) || total == INFINITY) return CRPX_ERROR_OVERFLOW;
  *value = total;
  return CRPX_OK;
}
static crpx_status
prepare (const crpx_clock_model *model, const crpx_tree *tree, const double *ages, crpx_tree_workspace *workspace,
         size_t count)
{
  if (!model || !tree || !ages || !workspace || !count || count != crpx_tree_node_count (tree) ||
      model->kind < CRPX_CLOCK_STRICT || model->kind > CRPX_CLOCK_AUTOCORRELATED || !isfinite (model->mean_rate) ||
      model->mean_rate <= 0 ||
      ((model->kind == CRPX_CLOCK_LOGNORMAL || model->kind == CRPX_CLOCK_AUTOCORRELATED) &&
       (!isfinite (model->dispersion) || model->dispersion <= 0)))
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_workspace_refresh (workspace, tree);
  if (status != CRPX_OK) return status;
  if (!crpx_tree_workspace_preorder (workspace)) {
    crpx_tree_workspace_reset (workspace);
    return CRPX_ERROR_INVALID_ARGUMENT;
  }
  for (size_t node = 0; node < count; node++) {
    crpx_node_id parent = crpx_tree_parent (tree, (crpx_node_id) node);
    if (!isfinite (ages[node]) ||
        (parent != CRPX_NODE_NONE && (!isfinite (ages[parent] - ages[node]) || ages[parent] < ages[node]))) {
      crpx_tree_workspace_reset (workspace);
      return CRPX_ERROR_INVALID_ARGUMENT;
    }
    if (parent != CRPX_NODE_NONE && !isfinite (variance (model, tree, ages, (crpx_node_id) node))) {
      crpx_tree_workspace_reset (workspace);
      return CRPX_ERROR_OVERFLOW;
    }
  }
  return CRPX_OK;
}
static double
variance (const crpx_clock_model *model, const crpx_tree *tree, const double *ages, crpx_node_id node)
{
  double v = model->dispersion * model->dispersion;
  if (model->kind != CRPX_CLOCK_AUTOCORRELATED) return model->kind == CRPX_CLOCK_LOGNORMAL ? v : 0;
  crpx_node_id parent = crpx_tree_parent (tree, node), ancestor = crpx_tree_parent (tree, parent);
  double interval = (ages[parent] - ages[node]) / 2;
  if (ancestor != CRPX_NODE_NONE) interval += (ages[ancestor] - ages[parent]) / 2;
  return v * interval;
}
