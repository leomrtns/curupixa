/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Derived analytic HKY process: biomcmc-recomb/src/likelihood.c and guenomu/lib/phylogeny.c. */
#include "substitution_model_internal.h"
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static unsigned group(size_t state) { return state == 0 || state == 2 ? 0u : 1u; }
static bool valid_mode(crpx_branch_mode mode)
{ return mode == CRPX_BRANCH_FIXED || mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL; }

crpx_status crpx_substitution_model_f81(size_t states, const double *frequencies, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!states || !frequencies) return CRPX_ERROR_INVALID_ARGUMENT;
  if (states > SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  double maximum = 0;
  for (size_t i = 0; i < states; ++i) {
    if (!isfinite(frequencies[i]) || frequencies[i] <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
    if (frequencies[i] > maximum) maximum = frequencies[i];
  }
  crpx_substitution_model *m = calloc(1, sizeof(*m));
  if (!m) return CRPX_ERROR_ALLOCATION;
  m->states = states; m->pi = malloc(states*sizeof(double));
  if (!m->pi) { free(m); return CRPX_ERROR_ALLOCATION; }
  long double total = 0, variation = 0, prefix = 0;
  for (size_t i = 0; i < states; ++i) total += (long double)frequencies[i]/maximum;
  for (size_t i = 0; i < states; ++i) {
    m->pi[i] = (double)(((long double)frequencies[i]/maximum)/total);
    if (!m->pi[i]) { crpx_substitution_model_free(m); return CRPX_ERROR_OVERFLOW; }
    variation += 2*(long double)m->pi[i]*prefix;
    prefix += m->pi[i];
  }
  m->base_decay = states == 1 ? 0 : (double)(1/variation);
  if (!isfinite(m->base_decay)) { crpx_substitution_model_free(m); return CRPX_ERROR_OVERFLOW; }
  *out = m;
  return CRPX_OK;
}
crpx_status crpx_substitution_model_poisson(size_t states, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!states) return CRPX_ERROR_INVALID_ARGUMENT;
  if (states > SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  double *weights = malloc(states*sizeof(double));
  if (!weights) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < states; ++i) weights[i] = 1;
  crpx_status status = crpx_substitution_model_f81(states, weights, out);
  free(weights); return status;
}
crpx_status crpx_substitution_model_hky(const double frequencies[4], double kappa, crpx_substitution_model **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!isfinite(kappa) || kappa <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_substitution_model *m;
  crpx_status status = crpx_substitution_model_f81(4, frequencies, &m);
  if (status != CRPX_OK) return status;
  m->hky = true; m->kappa = kappa;
  m->group_mass[0] = m->pi[0]+m->pi[2]; m->group_mass[1] = m->pi[1]+m->pi[3];
  long double scale = kappa > 1 ? kappa : 1;
  long double k = (long double)kappa/scale, v = 1/scale;
  long double d = 2*(k*((long double)m->pi[0]*m->pi[2]+(long double)m->pi[1]*m->pi[3])+
                    v*(long double)m->group_mass[0]*m->group_mass[1]);
  m->base_decay = (double)(v/d);
  m->transition_decay = (double)(k/d);
  for (unsigned g = 0; g < 2; ++g)
    m->group_decay[g] = (double)((k*m->group_mass[g]+v*m->group_mass[1-g])/d);
  if (!isfinite(m->transition_decay) || !isfinite(m->base_decay) || !isfinite(m->group_decay[0]) || !isfinite(m->group_decay[1])) {
    crpx_substitution_model_free(m); return CRPX_ERROR_OVERFLOW;
  }
  *out = m; return CRPX_OK;
}
void crpx_substitution_model_free(crpx_substitution_model *m) { if (m) { free(m->pi); free(m); } }
size_t crpx_substitution_model_states(const crpx_substitution_model *m) { return m ? m->states : 0; }
const double *crpx_substitution_model_frequencies(const crpx_substitution_model *m) { return m ? m->pi : NULL; }

/* Return 1 minus a decay factor accurately near zero. Overflow of a*t is the
 * well-defined long-edge limit rather than an arithmetic error. */
static double decay_loss(double rate, double t, crpx_branch_mode mode)
{
  double x = rate*t;
  if (isinf(x)) return 1;
  return mode == CRPX_BRANCH_FIXED ? -expm1(-x) : x/(1+x);
}
static double decay_factor(double rate, double t, crpx_branch_mode mode)
{
  double x = rate*t;
  if (mode == CRPX_BRANCH_FIXED) return exp(-x);
  return isinf(x) ? (1/rate)/t : 1/(1+x);
}
/* Positive algebra or a factored short-time series avoids cancellation when
 * kappa is tiny. Work in long double so log kernels retain sub-double probabilities. */
static long double loss_ratio(long double x) { return isinf(x) ? 1 : x/(1+x); }
static long double log_denominator(long double rate, double t)
{
  long double x = rate*(long double)t;
  return isinf(x) ? logl(rate)+logl(t) : log1pl(x);
}
static long double hky_entry(const crpx_substitution_model *m, size_t i, size_t j, double t,
                             crpx_branch_mode mode)
{
  if (!t) return i == j ? 1 : 0;
  unsigned group_id = group(i);
  long double r = m->group_mass[group_id], h = m->group_mass[1-group_id];
  long double b = m->base_decay, c = m->transition_decay, a = h*b+r*c;
  long double y = b*(long double)t, x = a*(long double)t, value;
  if (group(i) != group(j)) {
    value = mode == CRPX_BRANCH_FIXED ? -expm1l(-y) : loss_ratio(y);
    return (long double)m->pi[j]*value;
  }
  if (mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) {
    value = expl(logl(c)+logl(t)-log_denominator(a,t)-log_denominator(b,t))+loss_ratio(x)*loss_ratio(y);
  } else if (x <= 0.5L && y <= 0.5L) {
    long double d = (c-b)*(long double)t, power_y = y, power_x = x, polynomial = 1, factorial = 1;
    value = c*(long double)t;
    for (unsigned n = 2; n <= 64; ++n) {
      polynomial = x*polynomial+power_y;
      power_y *= y; power_x *= x; factorial *= n;
      long double coefficient = r >= 0.5L ? (power_x-h*power_y)/r : power_y+d*polynomial;
      long double term = coefficient/factorial;
      value += n%2 ? term : -term;
    }
  } else if (c >= b) {
    value = -expm1l(-y)+expl(-y)*(-expm1l(-r*(c-b)*(long double)t))/r;
  } else if (r < 0.5L) {
    value = -expm1l(-y)-expl(-x)*(-expm1l(-r*(b-c)*(long double)t))/r;
  } else {
    value = (-expm1l(-x)-h*(-expm1l(-y)))/r;
  }
  if (value < 0) value = 0;
  value *= m->pi[j];
  if (i == j) value += mode == CRPX_BRANCH_FIXED ? expl(-x) : expl(-log_denominator(a,t));
  return value;
}
static void hky_matrix(const crpx_substitution_model *m, double t, crpx_branch_mode mode, double matrix[16])
{
  for (size_t i = 0; i < 4; ++i) for (size_t j = 0; j < 4; ++j)
    matrix[i*4+j] = (double)hky_entry(m,i,j,t,mode);
}
void crpx_private_hky_log_matrix(const crpx_substitution_model *m, double t, crpx_branch_mode mode,
                                 double matrix[16])
{
  for (size_t i = 0; i < 4; ++i) for (size_t j = 0; j < 4; ++j) {
    long double value = hky_entry(m,i,j,t,mode);
    if (value > 0 || !t || i == j) matrix[i*4+j] = (double)logl(value);
    else if (m->base_decay*t > 1e-8 || m->group_decay[group(i)]*t > 1e-8 || !isfinite(value)) {
      matrix[i*4+j] = NAN; /* caller reports an unrepresentable analytic kernel */
    } else {
      /* Only an underflowed off-diagonal enters here. Its tiny-time leading
       * terms suffice: first-order direct change, second-order transversions.
       * This also supports platforms where long double has double's range. */
      double logt = log(t), logb = log(m->base_decay), logg;
      if (group(i) != group(j)) logg = logb+logt;
      else {
        double first = log(m->transition_decay)+logt;
        double second = log(m->group_mass[1-group(i)])+2*(logb+logt);
        if (mode == CRPX_BRANCH_FIXED) second -= log(2.0);
        double maximum = first > second ? first : second;
        logg = maximum == -INFINITY ? -INFINITY : maximum+log1p(exp(-fabs(first-second)));
      }
      matrix[i*4+j] = log(m->pi[j])+logg;
    }
  }
}
crpx_status crpx_substitution_model_transition(const crpx_substitution_model *m, double t,
                                              crpx_branch_mode mode, double *matrix, size_t capacity)
{
  if (!m || !matrix || !isfinite(t) || t < 0 || !valid_mode(mode)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (n > SIZE_MAX/n) return CRPX_ERROR_OVERFLOW;
  if (capacity < n*n) return CRPX_ERROR_OUT_OF_RANGE;
  if (m->hky) {
    double temporary[16]; hky_matrix(m,t,mode,temporary);
    for (size_t i = 0; i < 16; ++i) if (!isfinite(temporary[i])) return CRPX_ERROR_OVERFLOW;
    memcpy(matrix,temporary,sizeof(temporary));
  }
  else {
    double loss = decay_loss(m->base_decay,t,mode), factor = decay_factor(m->base_decay,t,mode);
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j)
      matrix[i*n+j] = m->pi[j]*loss+(i == j ? factor : 0);
  }
  return CRPX_OK;
}
crpx_status crpx_substitution_model_apply(const crpx_substitution_model *m, double t,
                                         crpx_branch_mode mode, const double *child, double *parent,
                                         size_t capacity)
{
  if (!m || !child || !parent || child == parent || !isfinite(t) || t < 0 || !valid_mode(mode))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (capacity < n) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t i = 0; i < n; ++i) if (!isfinite(child[i]) || child[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  if (m->hky) {
    double matrix[16], tmp[4]; hky_matrix(m,t,mode,matrix);
    for (size_t i = 0; i < 4; ++i) {
      long double sum = 0;
      for (size_t j = 0; j < 4; ++j) sum += (long double)matrix[i*4+j]*child[j];
      tmp[i] = (double)sum;
      if (!isfinite(tmp[i])) return CRPX_ERROR_OVERFLOW;
    }
    memcpy(parent,tmp,sizeof(tmp));
  } else {
    long double sum = 0;
    for (size_t i = 0; i < n; ++i) sum += (long double)m->pi[i]*child[i];
    double loss = decay_loss(m->base_decay,t,mode), factor = decay_factor(m->base_decay,t,mode);
    for (size_t i = 0; i < n; ++i)
      if (!isfinite((double)(loss*sum+(long double)factor*child[i]))) return CRPX_ERROR_OVERFLOW;
    for (size_t i = 0; i < n; ++i) parent[i] = (double)(loss*sum+(long double)factor*child[i]);
  }
  return CRPX_OK;
}
crpx_status crpx_substitution_model_generator(const crpx_substitution_model *m, double *matrix, size_t capacity)
{
  if (!m || !matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = m->states;
  if (n > SIZE_MAX/n) return CRPX_ERROR_OVERFLOW;
  if (capacity < n*n) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t i = 0; i < n; ++i) {
    long double sum = 0;
    for (size_t j = 0; j < n; ++j) if (i != j) {
      double rate = m->base_decay;
      if (m->hky && group(i) == group(j)) {
        rate = m->transition_decay;
      }
      matrix[i*n+j] = rate*m->pi[j]; sum += matrix[i*n+j];
    }
    matrix[i*n+i] = -(double)sum;
  }
  return CRPX_OK;
}

