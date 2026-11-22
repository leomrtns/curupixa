/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file kmedoids.c
    Deterministic PAM BUILD and the shared-accumulator FastPAM1 SWAP equations of Schubert and Rousseeuw (2021).
    These are generic dissimilarity calculations; no Euclidean coordinates or random initialisation are assumed. */
#include "kmedoids.h"
#include "clustering_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static crpx_status medoid_build (crpx_distance_generator_t generator, size_t component, crpx_kmedoids *result,
                                 bool *selected, double *nearest);
static crpx_status medoid_assign (crpx_distance_generator_t generator, size_t component, crpx_kmedoids *result,
                                  double *second, long double *cost);

void
crpx_kmedoids_free (crpx_kmedoids *result)
{
  if (result) { free (result->medoids); free (result->labels); free (result->distances); free (result); }
}

crpx_status
crpx_kmedoids_run (crpx_distance_generator_t generator, size_t component, size_t k,
                   const size_t *initial, size_t max_swaps, crpx_kmedoids **out)
{ // Cache each point's nearest two representatives; evaluate all removals together for each new candidate.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t size = crpx_distance_generator_samples (generator);
  if (!generator || !k || k > size || component >= crpx_distance_generator_components (generator))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof (size_t) || size > SIZE_MAX / sizeof (double) || k > SIZE_MAX / sizeof (long double))
    return CRPX_ERROR_OVERFLOW;
  crpx_kmedoids *result = calloc (1, sizeof (*result));
  bool *selected = calloc (size, sizeof (*selected));
  double *second = malloc (size * sizeof (*second));
  long double *removal = malloc (k * sizeof (*removal));
  crpx_status status = CRPX_OK;
  if (!result || !selected || !second || !removal) { status = CRPX_ERROR_ALLOCATION; goto done; }
  result->samples = size; result->clusters = k;
  result->medoids = malloc (k * sizeof (*result->medoids));
  result->labels = malloc (size * sizeof (*result->labels));
  result->distances = malloc (size * sizeof (*result->distances));
  if (!result->medoids || !result->labels || !result->distances) { status = CRPX_ERROR_ALLOCATION; goto done; }
  if (initial) {
    for (size_t i = 0; i < k; i++) {
      if (initial[i] >= size || selected[initial[i]]) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      result->medoids[i] = initial[i]; selected[initial[i]] = true;
    }
  } else {
    status = medoid_build (generator, component, result, selected, result->distances);
    if (status != CRPX_OK) goto done;
  }
  long double cost;
  status = medoid_assign (generator, component, result, second, &cost);
  if (status != CRPX_OK) goto done;
  result->converged = cost == 0 || k == size;
  while (!result->converged && result->swaps < max_swaps) {
    long double best_change = 0;
    size_t best_candidate = SIZE_MAX, best_slot = SIZE_MAX;
    for (size_t candidate = 0; candidate < size; candidate++) if (!selected[candidate]) {
      memset (removal, 0, k * sizeof (*removal));
      long double shared = 0;
      for (size_t i = 0; i < size; i++) {
        double distance;
        status = crpx_cluster_distance (generator, i, candidate, component, &distance);
        if (status != CRPX_OK) goto done;
        long double first_change = fminl (0, (long double)distance - result->distances[i]);
        shared += first_change;
        removal[result->labels[i]] += fminl (distance, second[i]) - result->distances[i] - first_change;
      }
      for (size_t slot = 0; slot < k; slot++) if (shared + removal[slot] < best_change) {
        best_change = shared + removal[slot]; best_candidate = candidate; best_slot = slot;
      }
    }
    /* Ignore changes below summation round-off, scaled to the supplied dissimilarities rather than an absolute unit. */
    if (best_candidate == SIZE_MAX || best_change >= -16 * LDBL_EPSILON * cost * size) {
      result->converged = true; break;
    }
    selected[result->medoids[best_slot]] = false;
    result->medoids[best_slot] = best_candidate; selected[best_candidate] = true;
    status = medoid_assign (generator, component, result, second, &cost);
    if (status != CRPX_OK) goto done;
    result->swaps++;
    if (cost == 0) result->converged = true;
  }
  if (!isfinite (cost) || cost > DBL_MAX) { status = CRPX_ERROR_OVERFLOW; goto done; }
  result->cost = (double)cost;
done:
  free (selected); free (second); free (removal);
  if (status != CRPX_OK) crpx_kmedoids_free (result);
  else *out = result;
  return status;
}

static crpx_status
medoid_build (crpx_distance_generator_t generator, size_t component, crpx_kmedoids *result,
               bool *selected, double *nearest)
{ // Greedy PAM BUILD chooses the largest reduction in total distance at each step.
  size_t size = result->samples;
  for (size_t i = 0; i < size; i++) nearest[i] = INFINITY;
  for (size_t slot = 0; slot < result->clusters; slot++) {
    size_t best = SIZE_MAX;
    long double best_cost = INFINITY;
    for (size_t candidate = 0; candidate < size; candidate++) if (!selected[candidate]) {
      long double cost = 0;
      for (size_t i = 0; i < size; i++) {
        double distance;
        crpx_status status = crpx_cluster_distance (generator, i, candidate, component, &distance);
        if (status != CRPX_OK) return status;
        cost += fmin (nearest[i], distance);
      }
      if (cost < best_cost) { best_cost = cost; best = candidate; }
    }
    if (best == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
    result->medoids[slot] = best; selected[best] = true;
    for (size_t i = 0; i < size; i++) {
      double distance;
      crpx_status status = crpx_cluster_distance (generator, i, best, component, &distance);
      if (status != CRPX_OK) return status;
      if (distance < nearest[i]) nearest[i] = distance;
    }
  }
  return CRPX_OK;
}

static crpx_status
medoid_assign (crpx_distance_generator_t generator, size_t component, crpx_kmedoids *result,
                double *second, long double *cost)
{
  *cost = 0;
  for (size_t i = 0; i < result->samples; i++) {
    double first = INFINITY; second[i] = INFINITY;
    size_t owner = SIZE_MAX;
    for (size_t slot = 0; slot < result->clusters; slot++) {
      double distance;
      crpx_status status = crpx_cluster_distance (generator, i, result->medoids[slot], component, &distance);
      if (status != CRPX_OK) return status;
      bool better = distance < first || (distance == first &&
                    (i == result->medoids[slot] ||
                     (i != result->medoids[owner] && result->medoids[slot] < result->medoids[owner])));
      if (better) { second[i] = first; first = distance; owner = slot; }
      else if (distance < second[i]) second[i] = distance;
    }
    result->labels[i] = owner; result->distances[i] = first; *cost += first;
  }
  return CRPX_OK;
}
