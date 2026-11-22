/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT
 * Derived analytic HKY process: biomcmc-recomb/src/likelihood.c and guenomu/lib/phylogeny.c. */
/* Construct finite-state substitution models used by likelihood and ancestral reconstruction. */
#include "substitution_model_internal.h"
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static crpx_status normalise_transition4 (double matrix[16]);
static crpx_status rates_transition (const crpx_substitution_model *, double, crpx_branch_mode, double[16]);
static unsigned
group (size_t state)
{
  return state == 0 || state == 2 ? 0u : 1u;
}
static bool
valid_mode (crpx_branch_mode mode)
{
  return mode == CRPX_BRANCH_FIXED || mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL;
}

crpx_status
crpx_substitution_model_f81 (size_t states, const double *frequencies, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!states || !frequencies) return CRPX_ERROR_INVALID_ARGUMENT;
  if (states > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double maximum = 0;
  for (size_t i = 0; i < states; ++i) {
    if (!isfinite (frequencies[i]) || frequencies[i] <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
    if (frequencies[i] > maximum) maximum = frequencies[i];
  }
  crpx_substitution_model *m = calloc (1, sizeof (*m));
  if (!m) return CRPX_ERROR_ALLOCATION;
  m->states = states;
  m->pi = malloc (states * sizeof (double));
  if (!m->pi) {
    free (m);
    return CRPX_ERROR_ALLOCATION;
  }
  long double total = 0, variation = 0, prefix = 0;
  for (size_t i = 0; i < states; ++i)
    total += (long double) frequencies[i] / maximum;
  for (size_t i = 0; i < states; ++i) {
    m->pi[i] = (double) (((long double) frequencies[i] / maximum) / total);
    if (!m->pi[i]) {
      crpx_substitution_model_free (m);
      return CRPX_ERROR_OVERFLOW;
    }
    variation += 2 * (long double) m->pi[i] * prefix;
    prefix += m->pi[i];
  }
  m->base_decay = states == 1 ? 0 : (double) (1 / variation);
  if (!isfinite (m->base_decay)) {
    crpx_substitution_model_free (m);
    return CRPX_ERROR_OVERFLOW;
  }
  *out = m;
  return CRPX_OK;
}
crpx_status
crpx_substitution_model_poisson (size_t states, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!states) return CRPX_ERROR_INVALID_ARGUMENT;
  if (states > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double *weights = malloc (states * sizeof (double));
  if (!weights) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < states; ++i)
    weights[i] = 1;
  crpx_status status = crpx_substitution_model_f81 (states, weights, out);
  free (weights);
  return status;
}
crpx_status
crpx_substitution_model_hky (const double frequencies[4], double kappa, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!isfinite (kappa) || kappa <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_substitution_model *m;
  crpx_status status = crpx_substitution_model_f81 (4, frequencies, &m);
  if (status != CRPX_OK) return status;
  m->hky = true;
  m->kappa = kappa;
  m->group_mass[0] = m->pi[0] + m->pi[2];
  m->group_mass[1] = m->pi[1] + m->pi[3];
  long double scale = kappa > 1 ? kappa : 1;
  long double k = (long double) kappa / scale, v = 1 / scale;
  long double d = 2 * (k * ((long double) m->pi[0] * m->pi[2] + (long double) m->pi[1] * m->pi[3]) +
                       v * (long double) m->group_mass[0] * m->group_mass[1]);
  m->base_decay = (double) (v / d);
  m->transition_decay = (double) (k / d);
  for (unsigned g = 0; g < 2; ++g)
    m->group_decay[g] = (double) ((k * m->group_mass[g] + v * m->group_mass[1 - g]) / d);
  if (!isfinite (m->transition_decay) || !isfinite (m->base_decay) || !isfinite (m->group_decay[0]) ||
      !isfinite (m->group_decay[1])) {
    crpx_substitution_model_free (m);
    return CRPX_ERROR_OVERFLOW;
  }
  *out = m;
  return CRPX_OK;
}
void
crpx_substitution_model_free (crpx_substitution_model *m)
{
  if (m) {
    free (m->pi);
    free (m->rates);
    free (m);
  }
}
size_t
crpx_substitution_model_states (const crpx_substitution_model *m)
{
  return m ? m->states : 0;
}
const double *
crpx_substitution_model_frequencies (const crpx_substitution_model *m)
{
  return m ? m->pi : NULL;
}

/* Return 1 minus a decay factor accurately near zero. Overflow of a*t is the well-defined long-edge limit rather than
   an arithmetic error. */
static double
decay_loss (double rate, double t, crpx_branch_mode mode)
{
  double x = rate * t;
  if (isinf (x)) return 1;
  return mode == CRPX_BRANCH_FIXED ? -expm1 (-x) : x / (1 + x);
}
static double
decay_factor (double rate, double t, crpx_branch_mode mode)
{
  double x = rate * t;
  if (mode == CRPX_BRANCH_FIXED) return exp (-x);
  return isinf (x) ? (1 / rate) / t : 1 / (1 + x);
}
/* Positive algebra or a factored short-time series avoids cancellation when kappa is tiny. Work in long double so log
   kernels retain sub-double probabilities. */
static long double
loss_ratio (long double x)
{
  return isinf (x) ? 1 : x / (1 + x);
}
static long double
log_denominator (long double rate, double t)
{
  long double x = rate * (long double) t;
  return isinf (x) ? logl (rate) + logl (t) : log1pl (x);
}
static long double
hky_entry (const crpx_substitution_model *m, size_t i, size_t j, double t, crpx_branch_mode mode)
{
  if (!t) return i == j ? 1 : 0;
  unsigned group_id = group (i);
  long double r = m->group_mass[group_id], h = m->group_mass[1 - group_id];
  long double b = m->base_decay, c = m->transition_decay, a = h * b + r * c;
  long double y = b * (long double) t, x = a * (long double) t, value;
  if (group (i) != group (j)) {
    value = mode == CRPX_BRANCH_FIXED ? -expm1l (-y) : loss_ratio (y);
    return (long double) m->pi[j] * value;
  }
  if (mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) {
    value =
      expl (logl (c) + logl (t) - log_denominator (a, t) - log_denominator (b, t)) + loss_ratio (x) * loss_ratio (y);
  } else if (x <= 0.5L && y <= 0.5L) {
    long double d = (c - b) * (long double) t, power_y = y, power_x = x, polynomial = 1, factorial = 1;
    value = c * (long double) t;
    for (unsigned n = 2; n <= 64; ++n) {
      polynomial = x * polynomial + power_y;
      power_y *= y;
      power_x *= x;
      factorial *= n;
      long double coefficient = r >= 0.5L ? (power_x - h * power_y) / r : power_y + d * polynomial;
      long double term = coefficient / factorial;
      value += n % 2 ? term : -term;
    }
  } else if (c >= b) {
    value = -expm1l (-y) + expl (-y) * (-expm1l (-r * (c - b) * (long double) t)) / r;
  } else if (r < 0.5L) {
    value = -expm1l (-y) - expl (-x) * (-expm1l (-r * (b - c) * (long double) t)) / r;
  } else {
    value = (-expm1l (-x) - h * (-expm1l (-y))) / r;
  }
  if (value < 0) value = 0;
  value *= m->pi[j];
  if (i == j) value += mode == CRPX_BRANCH_FIXED ? expl (-x) : expl (-log_denominator (a, t));
  return value;
}
static void
hky_matrix (const crpx_substitution_model *m, double t, crpx_branch_mode mode, double matrix[16])
{
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j)
      matrix[i * 4 + j] = (double) hky_entry (m, i, j, t, mode);
}
void
crpx_private_hky_log_matrix (const crpx_substitution_model *m, double t, crpx_branch_mode mode, double matrix[16])
{
  for (size_t i = 0; i < 4; ++i)
    for (size_t j = 0; j < 4; ++j) {
      long double value = hky_entry (m, i, j, t, mode);
      if (value > 0 || !t || i == j)
        matrix[i * 4 + j] = (double) logl (value);
      else if (m->base_decay * t > 1e-8 || m->group_decay[group (i)] * t > 1e-8 || !isfinite (value)) {
        matrix[i * 4 + j] = NAN; /* caller reports an unrepresentable analytic kernel */
      } else {
        /* Only an underflowed off-diagonal enters here. Its tiny-time leading terms suffice: first-order direct change,
   second-order transversions. This also supports platforms where long double has double's range. */
        double logt = log (t), logb = log (m->base_decay), logg;
        if (group (i) != group (j))
          logg = logb + logt;
        else {
          double first = log (m->transition_decay) + logt;
          double second = log (m->group_mass[1 - group (i)]) + 2 * (logb + logt);
          if (mode == CRPX_BRANCH_FIXED) second -= log (2.0);
          double maximum = first > second ? first : second;
          logg = maximum == -INFINITY ? -INFINITY : maximum + log1p (exp (-fabs (first - second)));
        }
        matrix[i * 4 + j] = log (m->pi[j]) + logg;
      }
    }
}
crpx_status
crpx_substitution_model_transition (const crpx_substitution_model *m, double t, crpx_branch_mode mode, double *matrix,
                                    size_t capacity)
{
  if (!m || !matrix || !isfinite (t) || t < 0 || !valid_mode (mode)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (n > SIZE_MAX / n) return CRPX_ERROR_OVERFLOW;
  if (capacity < n * n) return CRPX_ERROR_OUT_OF_RANGE;
  if (m->rates) {
    double temporary[16];
    crpx_status status = rates_transition (m, t, mode, temporary);
    if (status != CRPX_OK) return status;
    memcpy (matrix, temporary, sizeof (temporary));
  } else if (m->hky) {
    double temporary[16];
    hky_matrix (m, t, mode, temporary);
    for (size_t i = 0; i < 16; ++i)
      if (!isfinite (temporary[i])) return CRPX_ERROR_OVERFLOW;
    memcpy (matrix, temporary, sizeof (temporary));
  } else {
    double loss = decay_loss (m->base_decay, t, mode), factor = decay_factor (m->base_decay, t, mode);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        matrix[i * n + j] = m->pi[j] * loss + (i == j ? factor : 0);
  }
  return CRPX_OK;
}
crpx_status
crpx_substitution_model_apply (const crpx_substitution_model *m, double t, crpx_branch_mode mode, const double *child,
                               double *parent, size_t capacity)
{
  if (!m || !child || !parent || child == parent || !isfinite (t) || t < 0 || !valid_mode (mode))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (capacity < n) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t i = 0; i < n; ++i)
    if (!isfinite (child[i]) || child[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  if (m->rates) {
    double matrix[16], temporary[4];
    crpx_status status = rates_transition (m, t, mode, matrix);
    if (status != CRPX_OK) return status;
    for (size_t i = 0; i < 4; i++) {
      long double value = 0;
      for (size_t j = 0; j < 4; j++)
        value += (long double) matrix[i * 4 + j] * child[j];
      temporary[i] = (double) value;
      if (!isfinite (temporary[i])) return CRPX_ERROR_OVERFLOW;
    }
    memcpy (parent, temporary, sizeof (temporary));
  } else if (m->hky) {
    double matrix[16], tmp[4];
    hky_matrix (m, t, mode, matrix);
    for (size_t i = 0; i < 4; ++i) {
      long double sum = 0;
      for (size_t j = 0; j < 4; ++j)
        sum += (long double) matrix[i * 4 + j] * child[j];
      tmp[i] = (double) sum;
      if (!isfinite (tmp[i])) return CRPX_ERROR_OVERFLOW;
    }
    memcpy (parent, tmp, sizeof (tmp));
  } else {
    long double sum = 0;
    for (size_t i = 0; i < n; ++i)
      sum += (long double) m->pi[i] * child[i];
    double loss = decay_loss (m->base_decay, t, mode), factor = decay_factor (m->base_decay, t, mode);
    for (size_t i = 0; i < n; ++i)
      if (!isfinite ((double) (loss * sum + (long double) factor * child[i]))) return CRPX_ERROR_OVERFLOW;
    for (size_t i = 0; i < n; ++i)
      parent[i] = (double) (loss * sum + (long double) factor * child[i]);
  }
  return CRPX_OK;
}
crpx_status
crpx_substitution_model_generator (const crpx_substitution_model *m, double *matrix, size_t capacity)
{
  if (!m || !matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (n > SIZE_MAX / n) return CRPX_ERROR_OVERFLOW;
  if (capacity < n * n) return CRPX_ERROR_OUT_OF_RANGE;
  if (m->rates) {
    memcpy (matrix, m->rates, 16 * sizeof (double));
    return CRPX_OK;
  }
  for (size_t i = 0; i < n; ++i) {
    long double sum = 0;
    for (size_t j = 0; j < n; ++j)
      if (i != j) {
        double rate = m->base_decay;
        if (m->hky && group (i) == group (j)) {
          rate = m->transition_decay;
        }
        matrix[i * n + j] = rate * m->pi[j];
        sum += matrix[i * n + j];
      }
    matrix[i * n + i] = -(double) sum;
  }
  return CRPX_OK;
}

crpx_status
crpx_substitution_model_rates4 (const double generator[16], const double root_weights[4], crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!generator || !root_weights) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < 4; i++) {
    long double sum = 0, scale = 0;
    for (size_t j = 0; j < 4; j++) {
      double value = generator[4 * i + j];
      if (!isfinite (value) || (i != j && value < 0) || (i == j && value > 0)) return CRPX_ERROR_INVALID_ARGUMENT;
      sum += value;
      scale += fabs (value);
    }
    if (fabsl (sum) > 1e-12L * scale) return CRPX_ERROR_INVALID_ARGUMENT;
  }
  crpx_substitution_model *model;
  crpx_status status = crpx_substitution_model_f81 (4, root_weights, &model);
  if (status != CRPX_OK) return status;
  model->rates = malloc (16 * sizeof (double));
  if (!model->rates) {
    crpx_substitution_model_free (model);
    return CRPX_ERROR_ALLOCATION;
  }
  memcpy (model->rates, generator, 16 * sizeof (double));
  *out = model;
  return CRPX_OK;
}
crpx_status
crpx_substitution_model_binary_pair (const double rates[8], const double root_weights[4], crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!rates) return CRPX_ERROR_INVALID_ARGUMENT;
  static const size_t source[] = {0, 0, 1, 1, 2, 2, 3, 3}, target[] = {1, 2, 0, 3, 0, 3, 1, 2};
  double generator[16] = {0};
  for (size_t i = 0; i < 8; i++) {
    if (!isfinite (rates[i]) || rates[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
    generator[4 * source[i] + target[i]] = rates[i];
    generator[5 * source[i]] -= rates[i];
  }
  return crpx_substitution_model_rates4 (generator, root_weights, out);
}
static void
multiply4 (const double *first, const double *second, double *out)
{
  double temporary[16];
  for (size_t i = 0; i < 4; i++)
    for (size_t j = 0; j < 4; j++) {
      long double value = 0;
      for (size_t k = 0; k < 4; k++)
        value += (long double) first[4 * i + k] * second[4 * k + j];
      temporary[4 * i + j] = (double) value;
    }
  memcpy (out, temporary, sizeof (temporary));
}
static crpx_status
rates_transition (const crpx_substitution_model *model, double length, crpx_branch_mode mode, double out[16])
{ /* A small dense kernel extends trait models without slowing the analytic models. */
  if (mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) {
    long double augmented[4][8] = {{0}};
    for (size_t i = 0; i < 4; i++)
      for (size_t j = 0; j < 4; j++) {
        augmented[i][j] = (i == j) - (long double) length * model->rates[4 * i + j];
        augmented[i][j + 4] = i == j;
      }
    for (size_t column = 0; column < 4; column++) {
      size_t pivot = column;
      for (size_t i = column + 1; i < 4; i++)
        if (fabsl (augmented[i][column]) > fabsl (augmented[pivot][column])) pivot = i;
      if (augmented[pivot][column] == 0 || !isfinite (augmented[pivot][column])) return CRPX_ERROR_OVERFLOW;
      for (size_t j = 0; j < 8; j++) {
        long double swap = augmented[column][j];
        augmented[column][j] = augmented[pivot][j];
        augmented[pivot][j] = swap;
      }
      long double divisor = augmented[column][column];
      for (size_t j = 0; j < 8; j++)
        augmented[column][j] /= divisor;
      for (size_t i = 0; i < 4; i++)
        if (i != column) {
          long double factor = augmented[i][column];
          for (size_t j = 0; j < 8; j++)
            augmented[i][j] -= factor * augmented[column][j];
        }
    }
    for (size_t i = 0; i < 4; i++)
      for (size_t j = 0; j < 4; j++)
        out[4 * i + j] = (double) augmented[i][j + 4];
  } else {
    double rate = 0;
    for (size_t i = 0; i < 4; i++)
      if (-model->rates[5 * i] > rate) rate = -model->rates[5 * i];
    double mean = rate * length;
    if (!isfinite (mean)) return CRPX_ERROR_OVERFLOW;
    unsigned squares = 0;
    while (mean > .5) {
      mean *= .5;
      squares++;
    }
    double matrix[16], power[16] = {0};
    for (size_t i = 0; i < 4; i++)
      for (size_t j = 0; j < 4; j++) {
        matrix[4 * i + j] = (rate ? model->rates[4 * i + j] / rate : 0) + (i == j);
        power[4 * i + j] = i == j;
        out[4 * i + j] = (i == j) * exp (-mean);
      }
    double weight = exp (-mean);
    for (unsigned term = 1; term <= 64; term++) {
      multiply4 (power, matrix, power);
      weight *= mean / term;
      for (size_t i = 0; i < 16; i++)
        out[i] += weight * power[i];
      if (weight < 1e-19) break;
    }
    for (unsigned i = 0; i < squares; i++) {
      multiply4 (out, out, out);
      crpx_status status = normalise_transition4 (out);
      if (status != CRPX_OK) return status;
    }
  }
  return normalise_transition4 (out);
}
static crpx_status
normalise_transition4 (double matrix[16])
{
  for (size_t row = 0; row < 4; row++) {
    double total = 0;
    for (size_t column = 0; column < 4; column++) {
      double value = matrix[4 * row + column];
      if (!isfinite (value) || value < -1e-12) return CRPX_ERROR_OVERFLOW;
      matrix[4 * row + column] = fmax (0, value);
      total += matrix[4 * row + column];
    }
    if (!isfinite (total) || total <= 0) return CRPX_ERROR_OVERFLOW;
    for (size_t column = 0; column < 4; column++)
      matrix[4 * row + column] /= total;
  }
  return CRPX_OK;
}

void
crpx_private_rates_log_matrix (const crpx_substitution_model *model, double length, crpx_branch_mode mode,
                               double out[16])
{
  double matrix[16];
  crpx_status status = rates_transition (model, length, mode, matrix);
  for (size_t i = 0; i < 16; i++)
    out[i] = status == CRPX_OK ? log (matrix[i]) : NAN;
}

crpx_status
crpx_substitution_model_gtr (const double frequencies[4], const double exchangeabilities[6],
                             crpx_substitution_model **out)
{ /* Construct the reversible generator and remove its otherwise redundant overall rate scale. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!frequencies || !exchangeabilities) return CRPX_ERROR_INVALID_ARGUMENT;
  double pi[4], q[16] = {0}, total = 0, mean = 0;
  for (size_t i = 0; i < 4; i++) {
    if (!isfinite (frequencies[i]) || frequencies[i] <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
    total += frequencies[i];
  }
  if (!isfinite (total)) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < 4; i++) pi[i] = frequencies[i] / total;
  size_t k = 0;
  for (size_t i = 0; i < 4; i++) for (size_t j = i + 1; j < 4; j++, k++) {
    if (!isfinite (exchangeabilities[k]) || exchangeabilities[k] <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
    q[4 * i + j] = exchangeabilities[k] * pi[j];
    q[4 * j + i] = exchangeabilities[k] * pi[i];
  }
  for (size_t i = 0; i < 4; i++) for (size_t j = 0; j < 4; j++) if (i != j) {
    q[5 * i] -= q[4 * i + j];
    mean += pi[i] * q[4 * i + j];
  }
  if (!isfinite (mean) || mean <= 0) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < 16; i++) q[i] /= mean;
  return crpx_substitution_model_rates4 (q, pi, out);
}

bool
crpx_substitution_model_is_reversible (const crpx_substitution_model *model)
{ /* Detailed balance with the actual root distribution is necessary for root-invariant likelihoods. */
  if (!model) return false;
  if (!model->rates) return true;
  for (size_t i = 0; i < 4; i++) for (size_t j = i + 1; j < 4; j++) {
    double a = model->pi[i] * model->rates[4 * i + j], b = model->pi[j] * model->rates[4 * j + i];
    if (fabs (a - b) > 1e-12 * fmax (a, b)) return false;
  }
  return true;
}
