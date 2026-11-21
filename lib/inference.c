/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/* Evaluate conditional DTL-history probabilities and acceptance decisions; include sampling and tree-search APIs. */
#include "inference.h"
#include <math.h>

crpx_status
crpx_dtl_minimum_log_probability (const crpx_tree *gene, const crpx_dtl_species *species,
                                  const crpx_node_id *tip_map, const crpx_dtl_costs *costs,
                                  crpx_dtl_workspace *workspace, crpx_dtl_history_probability *result)
{ /* Normalise one selected minimum-cost history against all histories of the supplied trees. */
  crpx_dtl_history_probability probability;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_dtl_evaluate (gene, species, tip_map, costs, CRPX_DTL_LOG_PARTITION,
                                          workspace, &probability.log_partition);
  if (status != CRPX_OK) return status;
  status = crpx_dtl_evaluate (gene, species, tip_map, costs, CRPX_DTL_MINIMUM_COST,
                              workspace, &probability.minimum_cost);
  if (status != CRPX_OK) return status;
  probability.log_probability = -probability.minimum_cost - probability.log_partition;
  *result = probability;
  return CRPX_OK;
}

crpx_status
crpx_metropolis_accept (double log_target_change, double log_reverse_forward,
                        double log_uniform, bool *accepted)
{
  if (!accepted || isnan (log_target_change) || isnan (log_reverse_forward) ||
      isnan (log_uniform) || log_uniform >= 0) return CRPX_ERROR_INVALID_ARGUMENT;
  long double ratio = (long double) log_target_change + log_reverse_forward;
  if (isnan (ratio)) return CRPX_ERROR_INVALID_ARGUMENT;
  *accepted = ratio >= 0 || (long double) log_uniform < ratio;
  return CRPX_OK;
}
