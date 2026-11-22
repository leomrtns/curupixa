/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "clock.h"
#include <math.h>
crpx_status
crpx_clock_lengths (const crpx_tree *tree, const double *ages, const double *rates, double mean_rate, double *lengths,
                    size_t count)
{ /* Validate before writing, so a failed proposal leaves the previous lengths intact. */
  if (!tree || !ages || !lengths || !count || count != crpx_tree_node_count (tree) ||
      (!rates && (!isfinite (mean_rate) || mean_rate <= 0)) || lengths == ages || lengths == rates)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate (tree);
  if (status != CRPX_OK) return status;
  for (size_t node = 0; node < count; node++) {
    if (!isfinite (ages[node])) return CRPX_ERROR_INVALID_ARGUMENT;
    crpx_node_id parent = crpx_tree_parent (tree, (crpx_node_id) node);
    if (parent == CRPX_NODE_NONE) continue;
    double rate = rates ? rates[node] : mean_rate, duration = ages[parent] - ages[node];
    if (!isfinite (rate) || rate <= 0 || !isfinite (duration) || duration < 0) return CRPX_ERROR_INVALID_ARGUMENT;
    if (!isfinite (duration * rate)) return CRPX_ERROR_OVERFLOW;
  }
  for (size_t node = 0; node < count; node++) {
    crpx_node_id parent = crpx_tree_parent (tree, (crpx_node_id) node);
    lengths[node] = parent == CRPX_NODE_NONE ? 0 : (ages[parent] - ages[node]) * (rates ? rates[node] : mean_rate);
  }
  return CRPX_OK;
}
crpx_status
crpx_tree_ages (const crpx_tree *tree, const double *tip_ages, double tolerance, crpx_tree_workspace *workspace,
                double *ages, size_t count)
{ /* Reconstruct times from descendants, checking every independent path to a tip. */
  if (!tree || !workspace || !ages || !count || count != crpx_tree_node_count (tree) || !isfinite (tolerance) ||
      tolerance < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_workspace_refresh (workspace, tree);
  if (status != CRPX_OK) return status;
  const crpx_node_id *order = crpx_tree_workspace_postorder (workspace);
  if (!order) {
    crpx_tree_workspace_reset (workspace);
    return CRPX_ERROR_INVALID_ARGUMENT;
  }
  for (size_t i = 0; i < count && status == CRPX_OK; i++) {
    crpx_node_id node = order[i], child = crpx_tree_first_child (tree, node);
    if (child == CRPX_NODE_NONE) {
      ages[node] = tip_ages ? tip_ages[node] : 0;
      if (!isfinite (ages[node])) status = CRPX_ERROR_INVALID_ARGUMENT;
      continue;
    }
    bool first = true;
    for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling (tree, child)) {
      double length;
      status = crpx_tree_branch_length (tree, child, &length);
      if (status != CRPX_OK) break;
      double candidate = ages[child] + length;
      if (length < 0 || !isfinite (candidate)) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        break;
      }
      if (!first && fabs (candidate - ages[node]) > tolerance * fmax (1, fmax (fabs (candidate), fabs (ages[node])))) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        break;
      }
      ages[node] = candidate;
      first = false;
    }
  }
  crpx_tree_workspace_reset (workspace);
  return status;
}
