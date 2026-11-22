/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file stochastic_mapping.h
    \brief Sample CTMC branch histories conditional on endpoint states by uniformisation.
    Hobolth and Stone (2009), doi:10.1214/09-AOAS247. */
#ifndef CURUPIXA_STOCHASTIC_MAPPING_H
#define CURUPIXA_STOCHASTIC_MAPPING_H
#include "sampling.h"
#include "substitution_model.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_mapping_workspace crpx_mapping_workspace;
typedef struct {
  size_t count;
  double *times;
  size_t *states;
  double omitted_probability_bound;
} crpx_branch_history;
crpx_status crpx_mapping_workspace_new (crpx_mapping_workspace **out);
void crpx_mapping_workspace_free (crpx_mapping_workspace *workspace);
void crpx_branch_history_free (crpx_branch_history *history);
/* Times are sorted real state-change times in (0,length); states has count+1
   entries including the initial state. Self transitions used by uniformisation
   are removed. Conditional truncation error is bounded by tolerance; if the
   supplied virtual-event limit cannot achieve it, return CONVERGENCE, no history.
   Fixed branch lengths only. Sample ancestral endpoint states with
   crpx_likelihood_joint_states(sample=true), then call this for each branch.
   Numerical arrays are retained in workspace; result arrays own their memory. */
crpx_status crpx_substitution_bridge_sample (const crpx_substitution_model *model, double length, size_t first,
                                             size_t last, double tolerance, size_t max_virtual_events,
                                             crpx_uniform_function uniform, void *state,
                                             crpx_mapping_workspace *workspace, crpx_branch_history **out);
#ifdef __cplusplus
}
#endif
#endif
