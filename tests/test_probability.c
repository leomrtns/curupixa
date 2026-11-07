/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "prob_distribution.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; } } while (0)
#define CLOSE(a,b,tolerance) CHECK(fabs((a) - (b)) <= (tolerance))

int main(void)
{
  CLOSE(crpx_dnorm(0, 0, 1, false), 0.39894228040143267794, 1e-16);
  CLOSE(crpx_pnorm(0, 0, 1, false), 0.5, 1e-16);
  CLOSE(crpx_qnorm(0.975, 0, 1, false), 1.959963984540054, 1e-14);
  CLOSE(crpx_pnorm(-10, 0, 1, true), -53.23128515051247, 1e-12);
  CHECK(isnan(crpx_pnorm(NAN, 0, 1, false)));
  CHECK(isnan(crpx_qgamma(0.5, INFINITY, 1, false)));
  CHECK(isnan(crpx_qlnorm(0, 0, -1, false)));
  CHECK(isnan(crpx_qlnorm(1, 0, -1, false)));
  CHECK(isnan(crpx_qlnorm(-INFINITY, 0, -1, true)));
  CHECK(isnan(crpx_qlnorm(0, 0, -1, true)));
  CHECK(isnan(crpx_qnorm(2, 0, 0, false)));
  CHECK(isnan(crpx_qnorm(1, 0, 0, true)));
  CHECK(isnan(crpx_dpois(0, -1, false)));
  for (int i = -5; i <= 5; ++i) {
    double probability = crpx_pnorm(i, 1, 2, false);
    CLOSE(crpx_qnorm(probability, 1, 2, false), i, 1e-11);
    CLOSE(crpx_qnorm(log(probability), 1, 2, true), i, 1e-11);
  }
  double shapes[] = {0.1, 0.5, 1, 2, 10, 100};
  double probabilities[] = {1e-8, 0.001, 0.1, 0.5, 0.9, 0.999};
  for (size_t i = 0; i < sizeof(shapes) / sizeof(shapes[0]); ++i)
    for (size_t j = 0; j < sizeof(probabilities) / sizeof(probabilities[0]); ++j) {
      double quantile = crpx_qgamma(probabilities[j], shapes[i], 2, false);
      CLOSE(crpx_pgamma(quantile, shapes[i], 2, false), probabilities[j], 2e-12);
      CLOSE(crpx_pgamma(quantile, shapes[i], 2, true), log(probabilities[j]), 2e-9);
    }
  for (int i = 1; i < 20; ++i) {
    double x = i / 10.0;
    CLOSE(crpx_pgamma(x, 1, 2, false), -expm1(-2 * x), 2e-15);
    CLOSE(crpx_dgamma(x, 1, 2, false), 2 * exp(-2 * x), 2e-15);
    CLOSE(crpx_pgamma(x, 0.5, 1, false), erf(sqrt(x)), 2e-15);
  }
  CLOSE(crpx_ppois(0, 1000, true), -1000, 1e-12);
  CLOSE(crpx_ppois(0, 2, false), exp(-2), 1e-15);
  for (int lambda = 1; lambda <= 100; lambda += 11) {
    long double term = expl(-lambda), cumulative = 0;
    for (int x = 0; x < lambda + 20; ++x) {
      if (x) term *= (long double)lambda / x;
      cumulative += term;
      CLOSE(crpx_ppois(x, lambda, false), (double)cumulative, 1e-13);
      CLOSE(crpx_dpois(x, lambda, false), (double)term, 1e-14);
    }
    for (size_t j = 0; j < sizeof(probabilities) / sizeof(probabilities[0]); ++j) {
      double quantile = crpx_qpois(probabilities[j], lambda, false);
      CHECK(quantile >= 0 && quantile == floor(quantile));
      CHECK(crpx_ppois(quantile, lambda, false) >= probabilities[j] - 1e-14);
      CHECK(quantile == 0 || crpx_ppois(quantile - 1, lambda, false) < probabilities[j]);
      CHECK(crpx_qpois(log(probabilities[j]), lambda, true) == quantile);
    }
  }
  double rates[8], total = 0;
  CHECK(crpx_discrete_gamma(0.5, 0.5, rates, 8) == CRPX_OK);
  for (size_t i = 0; i < 8; ++i) { CHECK(rates[i] > 0); total += rates[i]; }
  CLOSE(total, 8, 1e-13);
  CHECK(crpx_discrete_gamma(1, DBL_MIN, rates, 8) == CRPX_OK);
  for (size_t i = 0; i < 8; ++i) CHECK(rates[i] > 0 && isfinite(rates[i]));
  CHECK(crpx_discrete_gamma(-1, 1, rates, 8) == CRPX_ERROR_INVALID_ARGUMENT);
  for (int i = 0; i <= 9; ++i) {
    double probability = crpx_pexp_dt(i, 0.3, 9, false);
    CHECK(crpx_qexp_dt(probability - 1e-15, 0.3, 9, false) == i);
  }
  CHECK(crpx_qexp_dt(1e-20, 0.3, 9, false) == 0);
  CHECK(crpx_qexp_dt(-1000, 0.3, 9, true) == 0);
  CLOSE(crpx_logspace_add(1000, 1000), 1000 + log(2), 1e-13);
  CHECK(crpx_logspace_add(-INFINITY, -INFINITY) == -INFINITY);
  CHECK(crpx_logspace_sub(2, 2) == -INFINITY);
  CHECK(isnan(crpx_logspace_sub(1, 2)));
  int sign;
  CLOSE(crpx_lgammafn(-0.5, &sign), log(2 * sqrt(3.14159265358979323846)), 1e-14);
  CHECK(sign == -1);

  crpx_rng_stream stream, duplicate;
  crpx_rng_stream_init(&stream, 11, 42);
  crpx_rng_stream_init(&duplicate, 11, 42);
  for (int i = 0; i < 100; ++i)
    CHECK(crpx_rng_gamma(&stream, 0.5, 2) == crpx_rng_gamma(&duplicate, 0.5, 2));
  double normal_sum = 0, gamma_sum = 0, poisson_sum = 0;
  for (int i = 0; i < 10000; ++i) {
    normal_sum += crpx_rng_norm(&stream, 2, 1);
    gamma_sum += crpx_rng_gamma(&stream, 0.5, 2);
    poisson_sum += crpx_rng_pois(&stream, 10);
  }
  CLOSE(normal_sum / 10000, 2, 0.05);
  CLOSE(gamma_sum / 10000, 0.25, 0.025);
  CLOSE(poisson_sum / 10000, 10, 0.15);

  double weights[] = {0, DBL_MAX, DBL_MAX, 0};
  crpx_discrete_sampler_t sampler = NULL;
  CHECK(crpx_discrete_sampler_new(weights, 4, &sampler) == CRPX_OK);
  CHECK(crpx_discrete_sampler_size(sampler) == 4);
  CHECK(crpx_discrete_sampler_pdf(sampler, 1, &total) == CRPX_OK && total == 0.5);
  size_t first = 0;
  for (int i = 0; i < 10000; ++i) {
    size_t index;
    CHECK(crpx_discrete_sampler_sample(sampler, &stream, &index) == CRPX_OK);
    CHECK(index == 1 || index == 2);
    first += index == 1;
  }
  CHECK(first > 4750 && first < 5250);
  crpx_discrete_sampler_free(sampler);
  double zeros[] = {0, 0};
  CHECK(crpx_discrete_sampler_new(zeros, 2, &sampler) == CRPX_ERROR_INVALID_ARGUMENT && !sampler);
  return 0;
}
