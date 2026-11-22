/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file phylogenetic_regression.h
    \brief Construct Brownian/OU trait covariances for the shared Gaussian regression API.
    Felsenstein (1985), doi:10.1086/284325; Pagel (1999), doi:10.1038/44766. */
#ifndef CURUPIXA_PHYLOGENETIC_REGRESSION_H
#define CURUPIXA_PHYLOGENETIC_REGRESSION_H
#include "tree.h"
#include "gaussian.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { CRPX_TRAIT_BROWNIAN, CRPX_TRAIT_OU_STATIONARY } crpx_trait_covariance_kind;
typedef struct {
  crpx_trait_covariance_kind kind;
  double alpha, lambda;
  const double *tip_error_variances;
} crpx_trait_covariance_options;
/* Default Brownian, lambda=1, no measurement error. OU uses a stationary root:
   Cij=exp(-alpha*path_distance)/(2alpha), alpha>0. Brownian Cij=root-to-MRCA length.
   lambda in [0,1] scales off-diagonals only. Optional errors add to the diagonal
   in units relative to the overall fitted diffusion variance. Tips are explicitly
   ordered, distinct leaves; subsets allow omission of missing observations without
   changing the tree. Errors and lengths must be finite and nonnegative.
   Reuse the returned factor in crpx_gaussian_regression() across markers/responses. */
void crpx_trait_covariance_options_init (crpx_trait_covariance_options *options);
crpx_status crpx_phylogenetic_covariance (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                                          const crpx_trait_covariance_options *options, crpx_covariance **out);
/* Prepare conditional ancestral moments for every node (including tips) using
   the same Brownian/OU covariance definition. Call crpx_gaussian_predict() with
   the chosen mean vectors. Brownian root mean is fixed by that vector; this is
   different from brownian_fit()'s flat-root, fitted-rate reconstruction. Observed
   error variances affect conditioning, but not the latent target variances. */
crpx_status crpx_phylogenetic_prediction (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                                          const crpx_trait_covariance_options *options, crpx_gaussian_prediction **out);
#ifdef __cplusplus
}
#endif
#endif
