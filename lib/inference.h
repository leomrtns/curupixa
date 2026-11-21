/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file inference.h
    \brief Evaluate conditional DTL-history probabilities and acceptance decisions; include sampling and tree-search APIs.
 */
#ifndef CURUPIXA_INFERENCE_H
#define CURUPIXA_INFERENCE_H
#include "dtl.h"
#include "amalgamation.h"
#include "tree_search.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double minimum_cost, log_partition, log_probability;
} crpx_dtl_history_probability;
/* Conditional probability of one selected minimum-cost history under the finite distribution P(history | observed gene
   tree, species tree, costs) proportional to exp(-history_cost). Reuses the DTL recurrence, first summing histories
   then minimising. log_probability=-minimum_cost-log_partition. Tied minimum histories each have this probability;
   their probabilities are not added together. This is not a gene-tree likelihood or a continuous-time birth-death
   model. The workspace finishes in minimum-cost mode, allowing explicit traceback. Outputs stay unchanged on error.
   Two dynamic-programming evaluations occur; allocated DP capacity is reused between them and across calls. */
crpx_status crpx_dtl_minimum_log_probability (const crpx_tree *gene, const crpx_dtl_species *species,
                                             const crpx_node_id *tip_map, const crpx_dtl_costs *costs,
                                             crpx_dtl_workspace *workspace,
                                             crpx_dtl_history_probability *result);
/* Caller supplies target log-density difference and log q(reverse)-log q(forward), including any Jacobian, and log(U)
   for U in [0,1). No RNG or chain state is hidden here. Signed infinities are supported; undefined infinity arithmetic
   and NaNs return INVALID_ARGUMENT without changing accepted. */
crpx_status crpx_metropolis_accept (double log_target_change, double log_reverse_forward,
                                   double log_uniform, bool *accepted);
#ifdef __cplusplus
}
#endif
#endif
