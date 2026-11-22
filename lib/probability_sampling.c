/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/* Explicit-stream distribution samplers and owned Walker alias tables. */
#include "prob_distribution.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
static double
uniform_positive (crpx_rng_stream *stream)
{
  double value;
  do {
    value = crpx_rng_stream_double (stream);
  } while (value == 0);
  return value;
}

double
crpx_rng_norm (crpx_rng_stream *stream, double mu, double sigma)
{
  if (!stream || !isfinite (mu) || !isfinite (sigma) || sigma < 0) return NAN;
  if (sigma == 0) return mu;
  double x, y, square;
  do {
    x = 2 * crpx_rng_stream_double (stream) - 1;
    y = 2 * crpx_rng_stream_double (stream) - 1;
    square = x * x + y * y;
  } while (square == 0 || square >= 1);
  return mu + sigma * x * sqrt (-2 * log (square) / square);
}

double
crpx_rng_lnorm (crpx_rng_stream *stream, double meanlog, double sdlog)
{
  return exp (crpx_rng_norm (stream, meanlog, sdlog));
}

/* Marsaglia and Tsang, A Simple Method for Generating Gamma Variables (2000). All proposal state is local; shape < 1
   uses the standard power transformation. */
double
crpx_rng_gamma (crpx_rng_stream *stream, double alpha, double beta)
{
  if (!stream || !isfinite (alpha) || !isfinite (beta) || alpha < 0 || beta <= 0) return NAN;
  if (alpha == 0) return 0;
  double shape = alpha < 1 ? alpha + 1 : alpha;
  double d = shape - 1.0 / 3.0, c = 1 / sqrt (9 * d);
  for (;;) {
    double normal = crpx_rng_norm (stream, 0, 1), v = 1 + c * normal;
    if (v <= 0) continue;
    v = v * v * v;
    double uniform = uniform_positive (stream), normal_squared = normal * normal;
    if (uniform < 1 - 0.0331 * normal_squared * normal_squared ||
        log (uniform) < 0.5 * normal_squared + d * (1 - v + log (v))) {
      double value = d * v / beta;
      return alpha < 1 ? value * pow (uniform_positive (stream), 1 / alpha) : value;
    }
  }
}

double
crpx_rng_pois (crpx_rng_stream *stream, double lambda)
{
  /* Doubles represent every integer only up to 2^53. */
  if (!stream || !isfinite (lambda) || lambda < 0 || lambda > 0x1p52) return NAN;
  if (lambda == 0) return 0;
  if (lambda < 10) {
    double threshold = exp (-lambda), product = 1;
    unsigned count = 0;
    do {
      product *= uniform_positive (stream);
      count++;
    } while (product > threshold);
    return count - 1;
  }
  /* Hoermann's transformed rejection (PTRS), avoiding repeated inverse CDFs. */
  double root = sqrt (lambda), b = .931 + 2.53 * root, a = -.059 + .02483 * b;
  double inverse = 1.1239 + 1.1328 / (b - 3.4), fast = .9277 - 3.6224 / (b - 2);
  for (;;) {
    double u = crpx_rng_stream_double (stream) - .5, v = uniform_positive (stream), us = .5 - fabs (u);
    if (us == 0) continue;
    double k = floor ((2 * a / us + b) * u + lambda + .43);
    if (k < 0) continue;
    if (us >= .07 && v <= fast) return k;
    if (us < .013 && v > us) continue;
    if (log (v * inverse / (a / (us * us) + b)) <= crpx_dpois (k, lambda, true)) return k;
  }
}

struct crpx_discrete_sampler_struct {
  size_t size;
  size_t *aliases;
  double *thresholds, *probabilities;
};

void
crpx_discrete_sampler_free (crpx_discrete_sampler_t sampler)
{
  if (sampler) {
    free (sampler->aliases);
    free (sampler->thresholds);
    free (sampler->probabilities);
    free (sampler);
  }
}

crpx_status
crpx_discrete_sampler_new (const double *weights, size_t size, crpx_discrete_sampler_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!weights || !size) return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof (double) || size > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  double maximum = 0;
  for (size_t i = 0; i < size; ++i) {
    if (!isfinite (weights[i]) || weights[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
    maximum = fmax (maximum, weights[i]);
  }
  if (maximum == 0) return CRPX_ERROR_INVALID_ARGUMENT;
  long double total = 0;
  for (size_t i = 0; i < size; ++i)
    total += (long double) weights[i] / maximum;
  crpx_discrete_sampler_t sampler = calloc (1, sizeof (*sampler));
  if (!sampler) return CRPX_ERROR_ALLOCATION;
  sampler->size = size;
  sampler->aliases = malloc (size * sizeof (size_t));
  sampler->thresholds = malloc (size * sizeof (double));
  sampler->probabilities = malloc (size * sizeof (double));
  double *scaled = malloc (size * sizeof (double));
  size_t *small = malloc (size * sizeof (size_t)), *large = malloc (size * sizeof (size_t));
  if (!sampler->aliases || !sampler->thresholds || !sampler->probabilities || !scaled || !small || !large) {
    free (scaled);
    free (small);
    free (large);
    crpx_discrete_sampler_free (sampler);
    return CRPX_ERROR_ALLOCATION;
  }
  size_t small_count = 0, large_count = 0;
  for (size_t i = 0; i < size; ++i) {
    sampler->probabilities[i] = (double) (((long double) weights[i] / maximum) / total);
    scaled[i] = (double) ((((long double) weights[i] / maximum) / total) * size);
    if (scaled[i] < 1)
      small[small_count++] = i;
    else
      large[large_count++] = i;
  }
  while (small_count && large_count) {
    size_t low = small[--small_count], high = large[--large_count];
    sampler->thresholds[low] = scaled[low];
    sampler->aliases[low] = high;
    scaled[high] = (scaled[high] - 1) + scaled[low];
    if (scaled[high] < 1)
      small[small_count++] = high;
    else
      large[large_count++] = high;
  }
  while (large_count) {
    size_t index = large[--large_count];
    sampler->thresholds[index] = 1;
    sampler->aliases[index] = index;
  }
  while (small_count) {
    size_t index = small[--small_count];
    sampler->thresholds[index] = 1;
    sampler->aliases[index] = index;
  }
  free (scaled);
  free (small);
  free (large);
  *out = sampler;
  return CRPX_OK;
}

size_t
crpx_discrete_sampler_size (const struct crpx_discrete_sampler_struct *sampler)
{
  return sampler ? sampler->size : 0;
}

crpx_status
crpx_discrete_sampler_sample (const struct crpx_discrete_sampler_struct *sampler, crpx_rng_stream *stream, size_t *out)
{
  if (!sampler || !stream || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t column;
  crpx_status status = crpx_rng_stream_range (stream, sampler->size, &column);
  if (status != CRPX_OK) return status;
  *out = crpx_rng_stream_double (stream) < sampler->thresholds[column] ? (size_t) column : sampler->aliases[column];
  return CRPX_OK;
}

crpx_status
crpx_discrete_sampler_pdf (const struct crpx_discrete_sampler_struct *sampler, size_t index, double *out)
{
  if (!sampler || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (index >= sampler->size) return CRPX_ERROR_OUT_OF_RANGE;
  *out = sampler->probabilities[index];
  return CRPX_OK;
}
