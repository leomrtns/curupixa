/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file tree_prior.h
    \brief Sample and evaluate tree priors with explicit conditioning and random state.
    Kingman (1982), doi:10.1016/0304-4149(82)90011-4;
    Gernhard (2008), doi:10.1016/j.jtbi.2008.04.005. */
#ifndef CURUPIXA_TREE_PRIOR_H
#define CURUPIXA_TREE_PRIOR_H
#include "sampling.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {
  CRPX_TREE_PRIOR_UNIFORM,
  CRPX_TREE_PRIOR_YULE,
  CRPX_TREE_PRIOR_BIRTH_DEATH,
  CRPX_TREE_PRIOR_COALESCENT
} crpx_tree_prior_kind;
typedef struct {
  crpx_tree_prior_kind kind;
  double birth_rate, death_rate, origin_age, coalescent_scale;
} crpx_tree_prior;
/* Uniform means labelled rooted binary topologies (PDA), without lengths.
   Yule/BD condition on n contemporary tips and fixed origin_age, with complete
   sampling. birth_rate>0, death_rate>=0; subcritical and critical BD are supported.
   They return the reconstructed extant tree, excluding extinct lineages and root
   stem. Coalescent uses constant pairwise waiting-time scale (e.g. 2Ne for diploids).
   All outputs are rooted binary; labels are copied. Independent calls share no state. */
crpx_status crpx_tree_prior_sample (const crpx_tree_prior *prior, size_t leaves, const char *const *labels,
                                    crpx_uniform_function uniform, void *state, crpx_tree **out);
/* Density of labelled, non-plane trees with respect to ordered branching ages;
   child order contributes no factor. Inputs for dated priors must be ultrametric,
   have positive edges and contemporary tips. Uniform ignores lengths. No fitting
   of parameters is performed. Outputs stay unchanged on failure. */
crpx_status crpx_tree_prior_log_density (const crpx_tree_prior *prior, const crpx_tree *tree, double *value);
/* Uniform labelled binary resolutions conditional on supplied clades. Existing
   lengths/labels are retained; new edges have length zero. Unary constraints are
   rejected by this inference interface (the lower constructor can preserve them).
   The density is a topology probability, not a density for new branch lengths. */
crpx_status crpx_tree_resolution_sample (const crpx_tree *constraint, crpx_uniform_function uniform, void *state,
                                         crpx_tree **out);
crpx_status crpx_tree_resolution_log_density (const crpx_tree *constraint, const crpx_node_id *constraint_tips,
                                              const crpx_tree *candidate, const crpx_node_id *candidate_tips,
                                              size_t taxa, double *value);
#ifdef __cplusplus
}
#endif
#endif
