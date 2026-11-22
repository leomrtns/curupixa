/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "kmedoids.h"
#include "global/random_number.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); return 1; } } while (0)
struct points { const double *values; bool fail; };
static crpx_status
distance (void *data, size_t first, size_t second, double *values, size_t components)
{
  (void)components;
  struct points *points = data;
  if (points->fail) return CRPX_ERROR_IO;
  *values = fabs (points->values[first] - points->values[second]);
  return CRPX_OK;
}
static double
cost (const double *points, size_t n, const size_t *medoids, size_t k)
{
  double sum = 0;
  for (size_t i = 0; i < n; i++) {
    double nearest = INFINITY;
    for (size_t j = 0; j < k; j++) nearest = fmin (nearest, fabs (points[i] - points[medoids[j]]));
    sum += nearest;
  }
  return sum;
}
int
main (void)
{
  crpx_rng_stream random;
  crpx_rng_stream_init (&random, 7152, 0);
  for (size_t repeat = 0; repeat < 100; repeat++) {
    double points[12];
    for (size_t i = 0; i < 12; i++) points[i] = (double)(crpx_rng_stream_next (&random) % 53);
    struct points data = {points, false};
    crpx_distance_generator_t generator = NULL;
    CHECK (crpx_distance_generator_new (12, 1, distance, &data, false, &generator) == CRPX_OK);
    size_t initial[3] = {0, 4, 8}, brute_medoids[3] = {0, 4, 8};
    crpx_kmedoids *result = NULL;
    CHECK (crpx_kmedoids_run (generator, 0, 3, initial, 100, &result) == CRPX_OK && result->converged);
    /* Independently enumerate every candidate/removal at every PAM step from the same initial medoids. */
    for (;;) {
      double current = cost (points, 12, brute_medoids, 3), best = current;
      size_t best_candidate = SIZE_MAX, best_slot = SIZE_MAX;
      for (size_t candidate = 0; candidate < 12; candidate++) {
        bool selected = false;
        for (size_t j = 0; j < 3; j++) if (brute_medoids[j] == candidate) selected = true;
        if (selected) continue;
        for (size_t slot = 0; slot < 3; slot++) {
          size_t old = brute_medoids[slot]; brute_medoids[slot] = candidate;
          double next = cost (points, 12, brute_medoids, 3); brute_medoids[slot] = old;
          if (next < best) { best = next; best_candidate = candidate; best_slot = slot; }
        }
      }
      if (best_candidate == SIZE_MAX) break;
      brute_medoids[best_slot] = best_candidate;
    }
    CHECK (result->cost == cost (points, 12, brute_medoids, 3));
    CHECK (result->cost == cost (points, 12, result->medoids, 3));
    for (size_t i = 0; i < 3; i++) CHECK (result->labels[result->medoids[i]] == i);
    crpx_kmedoids_free (result);
    CHECK (crpx_kmedoids_run (generator, 0, 1, NULL, 100, &result) == CRPX_OK && result->converged);
    for (size_t i = 0; i < 12; i++) CHECK (result->cost <= cost (points, 12, &i, 1));
    crpx_kmedoids_free (result);
    CHECK (crpx_kmedoids_run (generator, 0, 12, NULL, 100, &result) == CRPX_OK && result->cost == 0);
    crpx_kmedoids_free (result);
    data.fail = true;
    CHECK (crpx_kmedoids_run (generator, 0, 3, NULL, 100, &result) == CRPX_ERROR_IO && !result);
    crpx_distance_generator_free (generator);
  }
  return 0;
}
