/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
#include "graph.h"
#include "global/random_number.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); return 1; } } while (0)
struct matrix { size_t size; double *values; };
static crpx_status
weight (void *data, size_t first, size_t second, double *value)
{
  struct matrix *matrix = data;
  *value = matrix->values[first * matrix->size + second];
  return CRPX_OK;
}

int
main (void)
{
  crpx_rng_stream random;
  crpx_rng_stream_init (&random, 6728, 0);
  for (size_t repeat = 0; repeat < 100; repeat++) {
    enum { N = 8 };
    crpx_graph_edge edges[N * N];
    double weights[N * N], floyd[N * N];
    size_t count = 0;
    for (size_t i = 0; i < N; i++) for (size_t j = 0; j < N; j++) weights[i * N + j] = i == j ? 0 : INFINITY;
    bool directed = repeat % 2;
    for (size_t i = 0; i < N; i++) for (size_t j = directed ? 0 : i + 1; j < N; j++) {
      uint64_t draw = crpx_rng_stream_next (&random);
      if (i == j || draw % 3 == 0) continue;
      double cost = (double)(draw % 15);
      edges[count++] = (crpx_graph_edge){i, j, cost}; weights[i * N + j] = cost;
      if (!directed) weights[j * N + i] = cost;
    }
    crpx_graph *graph = NULL;
    CHECK (crpx_graph_new (N, edges, count, directed, &graph) == CRPX_OK);
    memcpy (floyd, weights, sizeof (floyd));
    for (size_t k = 0; k < N; k++) for (size_t i = 0; i < N; i++) for (size_t j = 0; j < N; j++)
      if (floyd[i * N + k] + floyd[k * N + j] < floyd[i * N + j]) floyd[i * N + j] = floyd[i * N + k] + floyd[k * N + j];
    for (size_t source = 0; source < N; source++) {
      double distances[N]; size_t previous[N];
      CHECK (crpx_dijkstra (graph, source, distances, previous, N) == CRPX_OK);
      for (size_t i = 0; i < N; i++) {
        CHECK (distances[i] == floyd[source * N + i]);
        if (i == source || isinf (distances[i])) CHECK (previous[i] == SIZE_MAX);
        else {
          double length = 0; size_t steps = 0, current = i;
          while (current != source && steps++ < N) {
            CHECK (previous[current] < N);
            length += weights[previous[current] * N + current]; current = previous[current];
          }
          CHECK (current == source && length == distances[i]);
        }
      }
    }
    crpx_spanning_forest *prim = NULL, *kruskal = NULL, *dense = NULL;
    if (!directed) {
      CHECK (crpx_prim (graph, &prim) == CRPX_OK && crpx_kruskal (graph, &kruskal) == CRPX_OK);
      struct matrix matrix = {N, weights};
      CHECK (crpx_prim_dense (N, weight, &matrix, &dense) == CRPX_OK);
      CHECK (prim->total == kruskal->total && prim->total == dense->total);
      CHECK (prim->components == kruskal->components && prim->components == dense->components);
      CHECK (prim->edge_count + prim->components == N);
      crpx_spanning_forest_free (prim); crpx_spanning_forest_free (kruskal); crpx_spanning_forest_free (dense);
    } else CHECK (crpx_prim (graph, &prim) == CRPX_ERROR_INVALID_ARGUMENT && !prim);
    crpx_graph_free (graph);
  }
  const crpx_graph_edge edges[] = {{0,1,5}, {0,1,-2}, {1,2,1}, {0,2,9}, {3,3,-5}};
  crpx_graph *graph = NULL;
  crpx_spanning_forest *forest = NULL;
  CHECK (crpx_graph_new (5, edges, SIZE_MAX, false, &graph) == CRPX_ERROR_OVERFLOW && !graph);
  CHECK (crpx_graph_new (5, edges, 5, false, &graph) == CRPX_OK);
  CHECK (crpx_prim (graph, &forest) == CRPX_OK && forest->total == -1 && forest->components == 3);
  crpx_spanning_forest_free (forest);
  CHECK (crpx_kruskal (graph, &forest) == CRPX_OK && forest->total == -1 && forest->components == 3);
  crpx_spanning_forest_free (forest);
  double distance[5] = {-1}; size_t previous[5] = {99};
  CHECK (crpx_dijkstra (graph, 0, distance, previous, 5) == CRPX_ERROR_INVALID_ARGUMENT && distance[0] == -1);
  crpx_graph_free (graph);
  const crpx_graph_edge huge[] = {{0,1,DBL_MAX},{1,2,DBL_MAX}};
  CHECK (crpx_graph_new (3, huge, 2, false, &graph) == CRPX_OK);
  CHECK (crpx_dijkstra (graph, 0, distance, previous, 5) == CRPX_ERROR_OVERFLOW && previous[0] == 99);
  CHECK (crpx_prim (graph, &forest) == CRPX_ERROR_OVERFLOW && !forest);
  crpx_graph_free (graph);
  CHECK (crpx_graph_new (0, NULL, 0, false, &graph) == CRPX_OK);
  CHECK (crpx_prim (graph, &forest) == CRPX_OK && !forest->edge_count && !forest->components);
  crpx_spanning_forest_free (forest); crpx_graph_free (graph);
  CHECK (crpx_prim_dense (0, NULL, NULL, &forest) == CRPX_OK && !forest->vertices);
  crpx_spanning_forest_free (forest);
  return 0;
}
