/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "particle_weights.h"
#include <math.h>
#include <float.h>
static crpx_status check_weights (const double *weights, size_t count);
static double conditional_ess (const double *weights, const double *scores, size_t count, double delta);
crpx_status
crpx_particle_normalise (const double *logs, size_t count, double *weights, double *log_sum, double *ess)
{
  if (!logs || !count || !weights || !log_sum || !ess) return CRPX_ERROR_INVALID_ARGUMENT;
  double maximum = -INFINITY;
  for (size_t i = 0; i < count; i++) {
    if (isnan (logs[i]) || logs[i] == INFINITY) return CRPX_ERROR_INVALID_ARGUMENT;
    if (logs[i] > maximum) maximum = logs[i];
  }
  if (!isfinite (maximum)) return CRPX_ERROR_NOT_FOUND;
  long double total = 0, squares = 0;
  for (size_t i = 0; i < count; i++)
    total += expl ((long double) logs[i] - maximum);
  for (size_t i = 0; i < count; i++) {
    weights[i] = (double) (expl ((long double) logs[i] - maximum) / total);
    squares += (long double) weights[i] * weights[i];
  }
  *log_sum = maximum + (double) logl (total);
  *ess = (double) (1 / squares);
  return CRPX_OK;
}
crpx_status
crpx_particle_temperature (const double *weights, const double *scores, size_t count, double beta, double fraction,
                           double *next_beta)
{ /* Conditional ESS controls the next increment even when current weights are unequal. */
  if (!scores || !next_beta || !isfinite (beta) || beta < 0 || beta >= 1 || !isfinite (fraction) || fraction <= 0 ||
      fraction >= 1)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = check_weights (weights, count);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < count; i++)
    if (!isfinite (scores[i])) return CRPX_ERROR_INVALID_ARGUMENT;
  double low = beta, high = 1, target = fraction * count;
  if (conditional_ess (weights, scores, count, 1 - beta) >= target) {
    *next_beta = 1;
    return CRPX_OK;
  }
  for (unsigned iteration = 0; iteration < 64; iteration++) {
    double middle = low + (high - low) / 2;
    if (middle == low || middle == high) break;
    if (conditional_ess (weights, scores, count, middle - beta) >= target)
      low = middle;
    else
      high = middle;
  }
  if (low == beta) return CRPX_ERROR_CONVERGENCE;
  *next_beta = low;
  return CRPX_OK;
}
crpx_status
crpx_particle_resample (const double *weights, size_t count, size_t offspring, bool stratified,
                        crpx_uniform_function uniform, void *state, size_t *indices)
{
  if (!offspring || !indices || !uniform) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = check_weights (weights, count);
  if (status != CRPX_OK) return status;
  double draw = 0;
  if (!stratified) {
    status = crpx_sample_uniform (uniform, state, &draw);
    if (status != CRPX_OK) return status;
  }
  size_t source = 0;
  long double cumulative = weights[0];
  for (size_t i = 0; i < offspring; i++) {
    if (stratified) {
      status = crpx_sample_uniform (uniform, state, &draw);
      if (status != CRPX_OK) return status;
    }
    long double point = ((long double) i + draw) / offspring;
    while (source + 1 < count && (weights[source] == 0 || point >= cumulative))
      cumulative += weights[++source];
    indices[i] = source;
  }
  return CRPX_OK;
}
static crpx_status
check_weights (const double *weights, size_t count)
{
  if (!weights || !count) return CRPX_ERROR_INVALID_ARGUMENT;
  long double sum = 0;
  for (size_t i = 0; i < count; i++) {
    if (!isfinite (weights[i]) || weights[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
    sum += weights[i];
  }
  return fabsl (sum - 1) <= 64 * DBL_EPSILON * count ? CRPX_OK : CRPX_ERROR_INVALID_ARGUMENT;
}
static double
conditional_ess (const double *weights, const double *scores, size_t count, double delta)
{
  double maximum = -INFINITY;
  for (size_t i = 0; i < count; i++)
    if (weights[i] > 0 && scores[i] > maximum) maximum = scores[i];
  long double first = 0, second = 0;
  for (size_t i = 0; i < count; i++) {
    if (weights[i] == 0) continue;
    long double value = expl ((long double) delta * ((long double) scores[i] - maximum));
    first += weights[i] * value;
    second += weights[i] * value * value;
  }
  return (double) (count * first * first / second);
}
