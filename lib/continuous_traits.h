/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_CONTINUOUS_TRAITS_H
#define CURUPIXA_CONTINUOUS_TRAITS_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_brownian_result crpx_brownian_result;
typedef struct crpx_pic_result crpx_pic_result;
typedef struct { bool reconstruct; } crpx_brownian_options;
/* Default fits the root mean and diffusion variance without allocating ancestral
 * arrays or performing the outside-subtree pass. Reconstruction is opt-in. */
void crpx_brownian_options_init(crpx_brownian_options *options);
/* Explicit tip_nodes[i] -> values[i] mapping. Every tree leaf must occur exactly
 * once; names and internal labels are irrelevant. Missing/nonfinite observations
 * are rejected. Both routines require >=2 tips and complete finite positive
 * non-root edge lengths; the root stem is ignored. Inputs are only borrowed for
 * the duration of the call. Results own their arrays, use original node IDs,
 * and remain valid after destroying the tree. No tree is changed or cached.
 * Construction failures set *out=NULL; independent calls are reentrant.
 * Both algorithms use O(nodes) time and memory, with no dense tip covariance. */
crpx_status crpx_brownian_fit(const crpx_tree *tree, const crpx_node_id *tip_nodes,
                             const double *values, size_t count, const crpx_brownian_options *options,
                             crpx_brownian_result **out);
/* Brownian ML: root mean fitted freely, sigma2 = residual quadratic form / tips.
 * Constant data have a boundary estimate and return INVALID_ARGUMENT.
 * Optional ancestral marginals use a flat prior on the root state and plug-in
 * sigma2: the root state is integrated rather than fixed at root_mean. Diffusion
 * rate and tree uncertainty are excluded. Polytomies and unary nodes are supported.
 * sigma2 has trait^2 / branch-length units; variances have trait^2 units.
 * The likelihood is the Gaussian density of observed tips at fitted parameters. */
void crpx_brownian_result_free(crpx_brownian_result *result);
size_t crpx_brownian_result_nodes(const crpx_brownian_result *result);
size_t crpx_brownian_result_tips(const crpx_brownian_result *result);
double crpx_brownian_result_root_mean(const crpx_brownian_result *result);
double crpx_brownian_result_sigma2(const crpx_brownian_result *result);
double crpx_brownian_result_log_likelihood(const crpx_brownian_result *result);
/* Borrowed arrays indexed by node ID, valid until result_free; absent in score-
 * only mode (NULL and count=0). NULL scalar results return NAN. */
const double *crpx_brownian_result_means(const crpx_brownian_result *result, size_t *count);
const double *crpx_brownian_result_variances(const crpx_brownian_result *result, size_t *count);

/* Standardised phylogenetically independent contrasts on a strictly binary
 * rooted tree. Contrast sign follows first-child minus second-child order.
 * Contrasts use trait / sqrt(branch-length) units. Returned variance is the
 * pre-standardisation evolutionary variance (adjusted branch-length units).
 * A nonbinary tree returns UNSUPPORTED; constant observations are permitted. */
crpx_status crpx_pic_calculate(const crpx_tree *tree, const crpx_node_id *tip_nodes,
                              const double *values, size_t count, crpx_pic_result **out);
void crpx_pic_result_free(crpx_pic_result *result);
size_t crpx_pic_result_count(const crpx_pic_result *result);
/* Parallel arrays in postorder contrast order; all views last until free. */
const crpx_node_id *crpx_pic_result_nodes(const crpx_pic_result *result);
const double *crpx_pic_result_contrasts(const crpx_pic_result *result);
const double *crpx_pic_result_variances(const crpx_pic_result *result);
#ifdef __cplusplus
}
#endif
#endif
