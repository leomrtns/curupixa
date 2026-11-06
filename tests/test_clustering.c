/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "clustering_goptics.h"
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; } } while (0)

typedef struct { const double *points; size_t calls; bool fail; } callback_state;
static crpx_status distances(void *data, size_t first, size_t second, double *out, size_t components)
{
  callback_state *state = data;
  state->calls++;
  if (state->fail) return CRPX_ERROR_IO;
  out[0] = fabs(state->points[first] - state->points[second]);
  if (components > 1) out[1] = out[0] * out[0];
  return CRPX_OK;
}

static void *worker(void *data)
{
  crpx_distance_generator_t generator = data;
  for (size_t i = 0; i < 7; ++i)
    for (size_t j = 0; j < 7; ++j) {
      double value;
      if (crpx_distance_generator_get(generator, i, j, 0, &value) != CRPX_OK) return (void *)1;
    }
  return NULL;
}

int main(void)
{
  const double points[] = {0, 0.1, 0.2, 10, 10.1, 10.2, 50};
  callback_state state = {points, 0, false};
  crpx_distance_generator_t generator;
  CHECK(crpx_distance_generator_new(7, 2, distances, &state, true, &generator) == CRPX_OK);
  pthread_t threads[4];
  for (size_t i = 0; i < 4; ++i) CHECK(pthread_create(threads + i, NULL, worker, generator) == 0);
  for (size_t i = 0; i < 4; ++i) {
    void *result;
    CHECK(pthread_join(threads[i], &result) == 0 && !result);
  }
  CHECK(state.calls == 21);
  double values[2] = {-1, -1};
  CHECK(crpx_distance_generator_get_all(generator, 0, 2, values, 2) == CRPX_OK);
  CHECK(values[0] == 0.2 && fabs(values[1] - 0.04) < 1e-16 && state.calls == 21);
  CHECK(crpx_distance_generator_get(generator, 7, 0, 0, values) == CRPX_ERROR_OUT_OF_RANGE);
  CHECK(crpx_distance_generator_get(generator, 0, 0, 2, values) == CRPX_ERROR_OUT_OF_RANGE);
  CHECK(crpx_distance_generator_reset(generator) == CRPX_OK);
  state.fail = true;
  CHECK(crpx_distance_generator_get(generator, 0, 1, 0, values) == CRPX_ERROR_IO);
  state.fail = false;
  CHECK(crpx_distance_generator_get(generator, 0, 1, 0, values) == CRPX_OK);
  CHECK(state.calls == 23);
  crpx_optics_t optics;
  CHECK(crpx_optics_run(generator, 0, 3, 0.25, &optics) == CRPX_OK);
  CHECK(crpx_optics_size(optics) == 7);
  const size_t *order = crpx_optics_order(optics);
  bool seen[7] = {false};
  for (size_t i = 0; i < 7; ++i) { CHECK(order[i] < 7 && !seen[order[i]]); seen[order[i]] = true; }
  size_t labels[7], clusters;
  CHECK(crpx_optics_extract_clusters(optics, 0.25, labels, 7, &clusters) == CRPX_OK);
  CHECK(clusters == 2);
  CHECK(labels[0] == labels[1] && labels[1] == labels[2]);
  CHECK(labels[3] == labels[4] && labels[4] == labels[5] && labels[0] != labels[3]);
  CHECK(labels[6] == SIZE_MAX);
  CHECK(isinf(crpx_optics_core_distances(optics)[6]));
  CHECK(crpx_optics_extract_clusters(optics, 0.3, labels, 7, &clusters) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_optics_free(optics);
  CHECK(crpx_optics_run(generator, 0, 1, 0, &optics) == CRPX_OK);
  CHECK(crpx_optics_extract_clusters(optics, 0, labels, 7, &clusters) == CRPX_OK && clusters == 7);
  crpx_optics_free(optics);
  CHECK(crpx_optics_run(generator, 0, 8, 1, &optics) == CRPX_ERROR_INVALID_ARGUMENT && !optics);
  crpx_distance_generator_free(generator);
  CHECK(crpx_distance_generator_new(7, 1, distances, &state, false, &generator) == CRPX_OK);
  state.calls = 0;
  CHECK(crpx_distance_generator_get(generator, 0, 1, 0, values) == CRPX_OK);
  CHECK(crpx_distance_generator_get(generator, 0, 1, 0, values) == CRPX_OK && state.calls == 2);
  CHECK(crpx_distance_generator_set_callback(generator, NULL, &state) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_distance_generator_free(generator);
  const double border_points[] = {0, 1, 2};
  callback_state border_state = {border_points, 0, false};
  CHECK(crpx_distance_generator_new(3, 1, distances, &border_state, true, &generator) == CRPX_OK);
  CHECK(crpx_optics_run(generator, 0, 3, 1, &optics) == CRPX_OK);
  CHECK(isinf(crpx_optics_reachability_distances(optics)[0]));
  CHECK(crpx_optics_extract_clusters(optics, 1, labels, 3, &clusters) == CRPX_OK && clusters == 1);
  CHECK(labels[0] == 0 && labels[1] == 0 && labels[2] == 0);
  crpx_optics_free(optics);
  crpx_distance_generator_free(generator);
  return 0;
}
