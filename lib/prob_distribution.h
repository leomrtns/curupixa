/* Adapted from biomcmc-lib and the R mathematical library.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins and the R Core Team
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_PROB_DISTRIBUTION_H
#define CURUPIXA_PROB_DISTRIBUTION_H
#include <stdbool.h>
#include <stddef.h>
#include "status.h"
#include "global/random_number.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Pure numerical functions return NAN for invalid parameters or failure. They do
 * not log, terminate, or share mutable caches. CDFs use the lower tail. log_p means
 * log output for density/CDF functions and log input for quantile functions.
 * Gamma uses shape alpha and RATE beta, so E[X] = alpha/beta. */
double crpx_dgamma(double x, double alpha, double beta, bool log_p);
double crpx_pgamma(double x, double alpha, double beta, bool log_p);
double crpx_qgamma(double p, double alpha, double beta, bool log_p);
double crpx_dnorm(double x, double mu, double sigma, bool log_p);
double crpx_pnorm(double x, double mu, double sigma, bool log_p);
double crpx_qnorm(double p, double mu, double sigma, bool log_p);
double crpx_dlnorm(double x, double meanlog, double sdlog, bool log_p);
double crpx_plnorm(double x, double meanlog, double sdlog, bool log_p);
double crpx_qlnorm(double p, double meanlog, double sdlog, bool log_p);
double crpx_dpois(double x, double lambda, bool log_p);
double crpx_ppois(double x, double lambda, bool log_p);
/* Quantiles and sampling restrict lambda to <= 2^52 so integer search progresses. */
double crpx_qpois(double p, double lambda, bool log_p);
/* Discrete exponential support is the integers 0..maximum. */
double crpx_dexp_dt(double x, double lambda, double maximum, bool log_p);
double crpx_pexp_dt(double x, double lambda, double maximum, bool log_p);
double crpx_qexp_dt(double p, double lambda, double maximum, bool log_p);
/* Equal-probability categories represented by their conditional means. */
crpx_status crpx_discrete_gamma(double alpha, double beta, double *rates, size_t count);

double crpx_lgammafn(double x, int *sign);
double crpx_gammafn(double x);
double crpx_log1p(double x);
double crpx_log1pmx(double x);
double crpx_expm1(double x);
double crpx_logspace_add(double log_x, double log_y);
double crpx_logspace_sub(double log_x, double log_y);

/* Streams belong to one worker at a time. Sampling never accesses a global RNG. */
double crpx_rng_norm(crpx_rng_stream *stream, double mu, double sigma);
double crpx_rng_lnorm(crpx_rng_stream *stream, double meanlog, double sdlog);
double crpx_rng_gamma(crpx_rng_stream *stream, double alpha, double beta);
double crpx_rng_pois(crpx_rng_stream *stream, double lambda);

typedef struct crpx_discrete_sampler_struct *crpx_discrete_sampler_t;
/* Finite non-negative weights are normalised; at least one must be positive.
 * Independent samplers own their arrays. Concurrent reads with separate streams
 * are safe. Construction takes O(n); each sample takes O(1). */
crpx_status crpx_discrete_sampler_new(const double *weights, size_t size, crpx_discrete_sampler_t *out);
void crpx_discrete_sampler_free(crpx_discrete_sampler_t sampler);
size_t crpx_discrete_sampler_size(const struct crpx_discrete_sampler_struct *sampler);
crpx_status crpx_discrete_sampler_sample(const struct crpx_discrete_sampler_struct *sampler,
                                       crpx_rng_stream *stream, size_t *out);
crpx_status crpx_discrete_sampler_pdf(const struct crpx_discrete_sampler_struct *sampler, size_t index, double *out);

#ifdef __cplusplus
}
#endif
#endif
