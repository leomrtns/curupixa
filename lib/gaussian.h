/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file gaussian.h
    \brief Reuse a positive-definite covariance factor across generalised least-squares fits.
    Breslow and Clayton (1993), doi:10.2307/2290687 (PQL). */
#ifndef CURUPIXA_GAUSSIAN_H
#define CURUPIXA_GAUSSIAN_H
#include "status.h"
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
  size_t size;
  double *lower;
  double log_determinant;
} crpx_covariance;
typedef struct crpx_regression_workspace crpx_regression_workspace;
typedef struct {
  size_t observations, predictors;
  double *coefficients, *coefficient_covariance;
  double residual_sum_squares, variance, log_likelihood;
} crpx_regression_result;
/* Factor a symmetric positive-definite row-major covariance; owns a copy.
   This numerical primitive is independent of how the covariance was constructed. */
crpx_status crpx_covariance_new (const double *matrix, size_t size, crpx_covariance **out);
void crpx_covariance_free (crpx_covariance *covariance);
crpx_status crpx_regression_workspace_new (crpx_regression_workspace **out);
void crpx_regression_workspace_free (crpx_regression_workspace *workspace);
/* X is observations*predictors row-major, including an intercept column if wanted.
   Missing observations must be removed explicitly before covariance preparation.
   Fits y~N(X beta,variance*C). reml=false uses ML variance RSS/n, true RSS/(n-p)
   and the restricted log likelihood. Returns coefficient covariance at fitted variance.
   Householder QR avoids forming X' C^-1 X. Singular designs return INVALID_ARGUMENT.
   Constant residuals give variance=0 and log_likelihood=+INFINITY (boundary fit).
   Reuse workspace for many responses/designs; each returned result owns its arrays. */
crpx_status crpx_gaussian_regression (const crpx_covariance *covariance, const double *design, size_t predictors,
                                      const double *response, bool reml, crpx_regression_workspace *workspace,
                                      crpx_regression_result **out);
void crpx_regression_result_free (crpx_regression_result *result);
/* Solve C*x=input using the retained factor; input/output may alias. */
crpx_status crpx_covariance_solve (const crpx_covariance *covariance, const double *input, double *output,
                                   size_t count);
typedef struct {
  size_t observations, predictors, iterations;
  double *coefficients, *latent, *probabilities;
  double log_penalised_likelihood;
} crpx_logistic_result;
/* Binary logistic mixed regression by penalised quasi-likelihood (PQL), with
   supplied positive-definite Gaussian random-effect covariance. A phylogenetic
   covariance gives phylogenetic logistic regression; its scale is a parameter,
   not estimated here. This approximation is not an exact marginal likelihood.
   X includes any intercept. Outputs remain NULL on non-convergence/separation.
   Dense O(n^3) working fits reuse the ordinary Gaussian regression kernel. */
crpx_status crpx_logistic_mixed_fit (const crpx_covariance *random_covariance, const double *design, size_t predictors,
                                     const uint8_t *response, double tolerance, size_t max_iterations,
                                     crpx_regression_workspace *workspace, crpx_logistic_result **out);
void crpx_logistic_result_free (crpx_logistic_result *result);
typedef struct {
  size_t observations, targets;
  double *gain, *variances;
} crpx_gaussian_prediction;
/* cross_covariance is targets*observations; target_variances is its unconditional
   diagonal. Returns read-only conditional gains and marginal variances. */
crpx_status crpx_gaussian_prediction_new (const crpx_covariance *observed, const double *cross_covariance,
                                          const double *target_variances, size_t targets,
                                          crpx_gaussian_prediction **out);
void crpx_gaussian_prediction_free (crpx_gaussian_prediction *prediction);
/* Allocation-free prediction with supplied mean vectors. Arrays must not overlap;
   covariance already includes any diffusion scale. No parameters are fitted. */
crpx_status crpx_gaussian_predict (const crpx_gaussian_prediction *prediction, const double *response,
                                   const double *observed_means, const double *target_means, double *means,
                                   double *variances);
#ifdef __cplusplus
}
#endif
#endif
