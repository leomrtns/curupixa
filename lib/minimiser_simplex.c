/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Nelder, J.A. and Mead, R. (1965), Computer Journal 7:308-313.
 * Replaces the incomplete GSL-derived draft with an independent implementation. */
#include "minimiser_simplex.h"
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct crpx_simplex {
  size_t n, best;
  bool initialised;
  crpx_objective objective;
  void *parameters;
  double *vertices, *values, *centroid, *reflection, *trial, *staged_vertices, *staged_values;
};

crpx_status crpx_simplex_new(size_t n, crpx_objective objective, void *parameters, crpx_simplex **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!n || !objective) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n == SIZE_MAX || n > SIZE_MAX/(n+1) || n*(n+1) > SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  crpx_simplex *s = calloc(1, sizeof(*s));
  if (!s) return CRPX_ERROR_ALLOCATION;
  s->n = n; s->objective = objective; s->parameters = parameters;
  s->vertices = malloc(n*(n+1)*sizeof(double));
  s->staged_vertices = malloc(n*(n+1)*sizeof(double));
  s->values = malloc((n+1)*sizeof(double)); s->staged_values = malloc((n+1)*sizeof(double));
  s->centroid = malloc(n*sizeof(double)); s->reflection = malloc(n*sizeof(double)); s->trial = malloc(n*sizeof(double));
  if (!s->vertices || !s->staged_vertices || !s->values || !s->staged_values ||
      !s->centroid || !s->reflection || !s->trial) { crpx_simplex_free(s); return CRPX_ERROR_ALLOCATION; }
  *out = s;
  return CRPX_OK;
}

void crpx_simplex_free(crpx_simplex *s)
{
  if (!s) return;
  free(s->vertices); free(s->staged_vertices); free(s->values); free(s->staged_values);
  free(s->centroid); free(s->reflection); free(s->trial); free(s);
}
static void refresh_best(crpx_simplex *s)
{
  s->best = 0;
  for (size_t i = 1; i <= s->n; ++i) if (s->values[i] < s->values[s->best]) s->best = i;
}
static crpx_status evaluate(crpx_simplex *s, const double *point, double *value)
{
  for (size_t i = 0; i < s->n; ++i) if (!isfinite(point[i])) return CRPX_ERROR_OVERFLOW;
  crpx_status status = s->objective(point, s->n, s->parameters, value);
  if (status != CRPX_OK) return status;
  return isfinite(*value) ? CRPX_OK : CRPX_ERROR_INVALID_ARGUMENT;
}
crpx_status crpx_simplex_initialise(crpx_simplex *s, const double *point, const double *steps)
{
  if (!s || !point || !steps) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t j = 0; j < s->n; ++j)
    if (!isfinite(point[j]) || !isfinite(steps[j]) || !steps[j] || point[j]+steps[j] == point[j])
      return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i <= s->n; ++i) {
    double *vertex = s->staged_vertices+i*s->n;
    memcpy(vertex, point, s->n*sizeof(double));
    if (i) vertex[i-1] += steps[i-1];
    crpx_status status = evaluate(s, vertex, s->staged_values+i);
    if (status != CRPX_OK) return status;
  }
  memcpy(s->vertices, s->staged_vertices, s->n*(s->n+1)*sizeof(double));
  memcpy(s->values, s->staged_values, (s->n+1)*sizeof(double));
  s->initialised = true; refresh_best(s);
  return CRPX_OK;
}
static void accept(crpx_simplex *s, size_t index, const double *point, double value)
{
  memcpy(s->vertices+index*s->n, point, s->n*sizeof(double));
  s->values[index] = value; refresh_best(s);
}
crpx_status crpx_simplex_iterate(crpx_simplex *s)
{
  if (!s || !s->initialised) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t worst = 0, second;
  for (size_t i = 1; i <= s->n; ++i) if (s->values[i] > s->values[worst]) worst = i;
  second = worst ? 0 : 1;
  for (size_t i = 0; i <= s->n; ++i) if (i != worst && s->values[i] > s->values[second]) second = i;
  for (size_t j = 0; j < s->n; ++j) {
    long double centre = 0;
    for (size_t i = 0; i <= s->n; ++i) if (i != worst) centre += (long double)s->vertices[i*s->n+j]/s->n;
    s->centroid[j] = (double)centre;
    s->reflection[j] = (double)(2*centre-s->vertices[worst*s->n+j]);
  }
  double reflected, trial_value;
  crpx_status status = evaluate(s, s->reflection, &reflected);
  if (status != CRPX_OK) return status;
  if (reflected < s->values[s->best]) {
    for (size_t j = 0; j < s->n; ++j) s->trial[j] = (double)(2*(long double)s->reflection[j]-s->centroid[j]);
    status = evaluate(s, s->trial, &trial_value);
    if (status != CRPX_OK) return status;
    if (trial_value < reflected) accept(s, worst, s->trial, trial_value);
    else accept(s, worst, s->reflection, reflected);
  } else if (reflected < s->values[second]) accept(s, worst, s->reflection, reflected);
  else {
    bool outside = reflected < s->values[worst];
    const double *endpoint = outside ? s->reflection : s->vertices+worst*s->n;
    for (size_t j = 0; j < s->n; ++j) s->trial[j] = 0.5*s->centroid[j]+0.5*endpoint[j];
    status = evaluate(s, s->trial, &trial_value);
    if (status != CRPX_OK) return status;
    if (outside ? trial_value <= reflected : trial_value < s->values[worst]) accept(s, worst, s->trial, trial_value);
    else {
      for (size_t i = 0; i <= s->n; ++i) {
        double *vertex = s->staged_vertices+i*s->n;
        if (i == s->best) { memcpy(vertex, s->vertices+i*s->n, s->n*sizeof(double)); s->staged_values[i] = s->values[i]; }
        else {
          for (size_t j = 0; j < s->n; ++j) vertex[j] = 0.5*s->vertices[i*s->n+j]+0.5*s->vertices[s->best*s->n+j];
          status = evaluate(s, vertex, s->staged_values+i);
          if (status != CRPX_OK) return status;
        }
      }
      memcpy(s->vertices, s->staged_vertices, s->n*(s->n+1)*sizeof(double));
      memcpy(s->values, s->staged_values, (s->n+1)*sizeof(double)); refresh_best(s);
    }
  }
  return CRPX_OK;
}
double crpx_simplex_size(const crpx_simplex *s)
{
  double maximum = 0;
  if (!s || !s->initialised) return NAN;
  for (size_t i = 0; i <= s->n; ++i) {
    double distance = 0;
    for (size_t j = 0; j < s->n; ++j) distance = hypot(distance, s->vertices[i*s->n+j]-s->vertices[s->best*s->n+j]);
    if (distance > maximum) maximum = distance;
  }
  return maximum;
}
crpx_status crpx_simplex_run(crpx_simplex *s, size_t maximum, double xtol, double ftol, size_t *iterations)
{
  if (!s || !s->initialised || !iterations || !isfinite(xtol) || !isfinite(ftol) || xtol < 0 || ftol < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  *iterations = 0;
  for (;;) {
    double spread = 0;
    for (size_t i = 0; i <= s->n; ++i) {
      double difference = fabs(s->values[i]-s->values[s->best]);
      if (difference > spread) spread = difference;
    }
    if (crpx_simplex_size(s) <= xtol && spread <= ftol) return CRPX_OK;
    if (*iterations == maximum) return CRPX_ERROR_CONVERGENCE;
    crpx_status status = crpx_simplex_iterate(s);
    if (status != CRPX_OK) return status;
    ++*iterations;
  }
}
const double *crpx_simplex_best_point(const crpx_simplex *s)
{ return s && s->initialised ? s->vertices+s->best*s->n : NULL; }
double crpx_simplex_best_value(const crpx_simplex *s) { return s && s->initialised ? s->values[s->best] : NAN; }
size_t crpx_simplex_dimensions(const crpx_simplex *s) { return s ? s->n : 0; }
