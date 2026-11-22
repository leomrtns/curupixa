/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "hdbscan.h"
#include "clustering_goptics.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); return 1; } } while (0)
struct points { const double *values; size_t size, calls, visits; bool fail; };
static crpx_status
distance (void *data, size_t first, size_t second, double *values, size_t components)
{
  (void)components;
  struct points *points = data;
  points->calls++;
  *values = fabs (points->values[first] - points->values[second]);
  return points->fail ? CRPX_ERROR_IO : CRPX_OK;
}
static crpx_status
neighbours (void *data, size_t query, double epsilon, size_t *ids, double *distances, size_t capacity, size_t *count)
{ // A sorted one-dimensional index locates a radius range with binary search, then visits only that range.
  struct points *points = data;
  if (points->fail) return CRPX_ERROR_IO;
  double centre = points->values[query];
  size_t left = 0, right = points->size;
  while (left < right) {
    size_t middle = left + (right - left) / 2;
    if (points->values[middle] < centre - epsilon) left = middle + 1; else right = middle;
  }
  *count = 0;
  for (size_t i = left; i < points->size && points->values[i] <= centre + epsilon; i++) if (i != query) {
    if (*count >= capacity) return CRPX_ERROR_OUT_OF_RANGE;
    ids[*count] = i; distances[(*count)++] = fabs (points->values[i] - centre); points->visits++;
  }
  return CRPX_OK;
}
int
main (void)
{
  const double values[] = {0, 0.1, 0.2, 10, 10.1, 10.2, 50};
  struct points points = {values, 7, 0, 0, false};
  crpx_distance_generator_t generator = NULL;
  CHECK (crpx_distance_generator_new (7, 1, distance, &points, false, &generator) == CRPX_OK);
  crpx_optics_t scan = NULL, indexed = NULL;
  CHECK (crpx_optics_run (generator, 0, 3, 0.3, &scan) == CRPX_OK);
  CHECK (crpx_optics_run_neighbours (7, neighbours, &points, 3, 0.3, &indexed) == CRPX_OK);
  CHECK (points.calls == 42 && points.visits == 12);
  for (size_t i = 0; i < 7; i++) {
    CHECK (crpx_optics_order (scan)[i] == crpx_optics_order (indexed)[i]);
    CHECK (crpx_optics_core_distances (scan)[i] == crpx_optics_core_distances (indexed)[i]);
    CHECK (crpx_optics_reachability_distances (scan)[i] == crpx_optics_reachability_distances (indexed)[i]);
    CHECK (crpx_optics_predecessors (scan)[i] == crpx_optics_predecessors (indexed)[i]);
  }
  crpx_optics_free (scan); crpx_optics_free (indexed);
  crpx_hdbscan *result = NULL;
  CHECK (crpx_hdbscan_run (generator, 0, 2, 3, CRPX_HDBSCAN_EOM, false, &result) == CRPX_OK);
  CHECK (result->clusters == 2 && result->labels[6] == SIZE_MAX);
  CHECK (result->labels[0] == result->labels[1] && result->labels[1] == result->labels[2]);
  CHECK (result->labels[3] == result->labels[4] && result->labels[4] == result->labels[5]);
  for (size_t i = 0; i < 7; i++) CHECK (isfinite (result->probabilities[i]) && result->probabilities[i] >= 0 && result->probabilities[i] <= 1);
  crpx_hdbscan_free (result);
  points.fail = true;
  CHECK (crpx_hdbscan_run (generator, 0, 2, 3, CRPX_HDBSCAN_EOM, false, &result) == CRPX_ERROR_IO && !result);
  CHECK (crpx_optics_run_neighbours (7, neighbours, &points, 3, 0.3, &indexed) == CRPX_ERROR_IO && !indexed);
  crpx_distance_generator_free (generator);
  const double duplicates[] = {0, 0, 0, 10, 10, 10};
  points = (struct points){duplicates, 6, 0, 0, false};
  CHECK (crpx_distance_generator_new (6, 1, distance, &points, false, &generator) == CRPX_OK);
  CHECK (crpx_hdbscan_run (generator, 0, 2, 3, CRPX_HDBSCAN_EOM, false, &result) == CRPX_OK && result->clusters == 2);
  for (size_t i = 0; i < 6; i++) CHECK (result->labels[i] != SIZE_MAX && result->probabilities[i] == 1);
  crpx_hdbscan_free (result);
  const double identical[6] = {0}; points.values = identical;
  CHECK (crpx_hdbscan_run (generator, 0, 2, 2, CRPX_HDBSCAN_EOM, false, &result) == CRPX_OK && !result->clusters);
  crpx_hdbscan_free (result);
  CHECK (crpx_hdbscan_run (generator, 0, 2, 2, CRPX_HDBSCAN_EOM, true, &result) == CRPX_OK && result->clusters == 1);
  for (size_t i = 0; i < 6; i++) CHECK (result->labels[i] == 0 && result->probabilities[i] == 1);
  crpx_hdbscan_free (result);
  CHECK (crpx_hdbscan_run (generator, 0, 2, 2, CRPX_HDBSCAN_LEAF, true, &result) == CRPX_OK && !result->clusters);
  crpx_hdbscan_free (result); crpx_distance_generator_free (generator);
  crpx_graph_edge tied_edges[5];
  for (size_t i = 0; i < 5; i++) tied_edges[i] = (crpx_graph_edge){i, i + 1, 1};
  crpx_spanning_forest tied_mst = {6, 5, 1, tied_edges, 5};
  CHECK (crpx_hdbscan_from_mst (&tied_mst, 2, CRPX_HDBSCAN_EOM, false, &result) == CRPX_OK && !result->clusters);
  crpx_hdbscan_free (result);
  CHECK (crpx_hdbscan_from_mst (&tied_mst, 2, CRPX_HDBSCAN_EOM, true, &result) == CRPX_OK && result->clusters == 1);
  for (size_t i = 0; i < 6; i++) CHECK (result->probabilities[i] == 1);
  crpx_hdbscan_free (result);
  /* A large unbranched hierarchy exercises iterative condensation without a recursion-depth limit. */
  size_t size = 10000;
  crpx_graph_edge *edges = malloc ((size - 1) * sizeof (*edges));
  CHECK (edges);
  for (size_t i = 0; i < size - 1; i++) edges[i] = (crpx_graph_edge){i, i + 1, (double)(i + 1)};
  crpx_spanning_forest mst = {size, size - 1, 1, edges, 0};
  CHECK (crpx_hdbscan_from_mst (&mst, 3, CRPX_HDBSCAN_EOM, true, &result) == CRPX_OK);
  CHECK (result->clusters == 1 && result->edge_count == size);
  crpx_hdbscan_free (result);
  edges[0].first = edges[0].second;
  CHECK (crpx_hdbscan_from_mst (&mst, 3, CRPX_HDBSCAN_EOM, true, &result) == CRPX_ERROR_TOPOLOGY && !result);
  free (edges);
  return 0;
}
