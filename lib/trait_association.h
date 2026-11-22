/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file trait_association.h
    \brief Reusable binary-trait pairing and phylogenetic signal calculations.
    Brynildsrud et al. (2016), doi:10.1186/s13059-016-1108-8; Fritz and Purvis (2010), doi:10.1111/j.1523-1739.2010.01455.x. */
#ifndef CURUPIXA_TRAIT_ASSOCIATION_H
#define CURUPIXA_TRAIT_ASSOCIATION_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_association_workspace crpx_association_workspace;
typedef struct {
  size_t pairs, supporting_minimum, supporting_maximum, observed_tips;
} crpx_pairing_result;
crpx_status crpx_association_workspace_new (crpx_association_workspace **out);
void crpx_association_workspace_free (crpx_association_workspace *workspace);
/* Maximise edge-disjoint tip pairs differing in both binary characters; report
   support extremes among ALL maximum pairings. 00/11 supports positive association,
   01/10 opposes it. Strictly binary rooted tree; branch lengths ignored.
   tips lists every tree leaf exactly once; values are 0/1, UINT8_MAX excludes a tip.
   O(nodes) time and retained worker memory; no permutations or multiple-testing policy. */
crpx_status crpx_trait_contrasting_pairs (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                                          const uint8_t *first, const uint8_t *second,
                                          crpx_association_workspace *workspace, crpx_pairing_result *result);
/* Raw sum of absolute child-to-parent changes using inverse adjusted-branch-length
   nodal means, as in caper phylo.d. Requires complete binary observations and
   strictly positive lengths. Reuses the Brownian Gaussian-message primitive.
   This raw dispersion is NOT normalised D. */
crpx_status crpx_trait_binary_dispersion (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                                          const uint8_t *values, crpx_association_workspace *workspace, double *score);
/* D=(observed-brownian_mean)/(random_mean-brownian_mean). Reference means are
   explicit: a downstream program can obtain them by repeated calls above with
   shuffled tips and Brownian-threshold simulations. No permutation is mandatory
   inside the score function, but reference calibration is necessary to call it D. */
crpx_status crpx_trait_d_statistic (double observed, double brownian_mean, double random_mean, double *value);
#ifdef __cplusplus
}
#endif
#endif
