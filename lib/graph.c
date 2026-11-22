/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file graph.c
    Dijkstra shortest paths and Prim/Kruskal spanning forests, following the original graph algorithms.
    Dense Prim also supplies the mutual-reachability spanning tree used by HDBSCAN. */
#include "graph.h"
#include "priority_queue_internal.h"
#include "union_find_internal.h"
#include <float.h>
#include <math.h>
#include <string.h>

static crpx_status graph_check (const crpx_graph *graph, bool nonnegative);
static crpx_status forest_new (size_t vertices, crpx_spanning_forest **out);
static crpx_status forest_total (crpx_spanning_forest *forest);
static int edge_compare (const void *first, const void *second);

crpx_status
crpx_graph_new (size_t vertices, const crpx_graph_edge *edges, size_t count, bool directed, crpx_graph **out)
{ // Build compressed adjacency arrays once, so repeated searches need no edge-list conversion.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if ((!edges && count) || (!vertices && count)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (vertices == SIZE_MAX || vertices + 1 > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  if (count > SIZE_MAX / sizeof (crpx_graph_edge)) return CRPX_ERROR_OVERFLOW;
  size_t arcs = count;
  for (size_t i = 0; i < count; i++) {
    if (edges[i].first >= vertices || edges[i].second >= vertices || !isfinite (edges[i].weight))
      return CRPX_ERROR_INVALID_ARGUMENT;
    if (!directed && edges[i].first != edges[i].second) {
      if (arcs == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
      arcs++;
    }
  }
  if (arcs > SIZE_MAX / sizeof (double) || arcs > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  crpx_graph *graph = calloc (1, sizeof (*graph));
  size_t *cursor = vertices ? malloc (vertices * sizeof (*cursor)) : NULL;
  if (!graph || (vertices && !cursor)) { free (graph); free (cursor); return CRPX_ERROR_ALLOCATION; }
  graph->vertices = vertices; graph->arcs = arcs; graph->directed = directed;
  graph->offsets = calloc (vertices + 1, sizeof (*graph->offsets));
  graph->targets = arcs ? malloc (arcs * sizeof (*graph->targets)) : NULL;
  graph->weights = arcs ? malloc (arcs * sizeof (*graph->weights)) : NULL;
  if (!graph->offsets || (arcs && (!graph->targets || !graph->weights))) {
    crpx_graph_free (graph); free (cursor); return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < count; i++) {
    graph->offsets[edges[i].first + 1]++;
    if (!directed && edges[i].first != edges[i].second) graph->offsets[edges[i].second + 1]++;
  }
  for (size_t i = 0; i < vertices; i++) graph->offsets[i + 1] += graph->offsets[i];
  if (vertices) memcpy (cursor, graph->offsets, vertices * sizeof (*cursor));
  for (size_t i = 0; i < count; i++) {
    size_t a = edges[i].first, b = edges[i].second, position = cursor[a]++;
    graph->targets[position] = b; graph->weights[position] = edges[i].weight;
    if (!directed && a != b) {
      position = cursor[b]++; graph->targets[position] = a; graph->weights[position] = edges[i].weight;
    }
  }
  free (cursor); *out = graph;
  return CRPX_OK;
}

void
crpx_graph_free (crpx_graph *graph)
{
  if (graph) { free (graph->offsets); free (graph->targets); free (graph->weights); free (graph); }
}

void
crpx_spanning_forest_free (crpx_spanning_forest *forest)
{
  if (forest) { free (forest->edges); free (forest); }
}

crpx_status
crpx_dijkstra (const crpx_graph *graph, size_t source, double *distances, size_t *predecessors, size_t capacity)
{ // Finalise the nearest unsettled vertex, then relax its outgoing arcs.
  crpx_status status = graph_check (graph, true);
  if (status != CRPX_OK) return status;
  size_t size = graph->vertices;
  if (source >= size || !distances || !predecessors) return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity < size) return CRPX_ERROR_OUT_OF_RANGE;
  crpx_priority_queue queue;
  status = crpx_queue_init (&queue, size);
  if (status != CRPX_OK) return status;
  long double *lengths = malloc (size * sizeof (*lengths));
  size_t *previous = malloc (size * sizeof (*previous));
  bool *settled = calloc (size, sizeof (*settled));
  if (!lengths || !previous || !settled) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i = 0; i < size; i++) { lengths[i] = INFINITY; previous[i] = SIZE_MAX; }
  lengths[source] = 0; crpx_queue_offer (&queue, source, 0);
  while (queue.size) {
    size_t current = crpx_queue_pop (&queue);
    settled[current] = true;
    for (size_t i = graph->offsets[current]; i < graph->offsets[current + 1]; i++) {
      size_t next = graph->targets[i];
      if (settled[next]) continue;
      long double candidate = lengths[current] + graph->weights[i];
      if (candidate < lengths[next] || (candidate == INFINITY && previous[next] == SIZE_MAX)) {
        lengths[next] = candidate; previous[next] = current;
        crpx_queue_offer (&queue, next, candidate);
      }
    }
  }
  for (size_t i = 0; i < size; i++) if (previous[i] != SIZE_MAX && lengths[i] > DBL_MAX) {
    status = CRPX_ERROR_OVERFLOW; goto done;
  }
  for (size_t i = 0; i < size; i++) { distances[i] = (double)lengths[i]; predecessors[i] = previous[i]; }
done:
  free (lengths); free (previous); free (settled); crpx_queue_free (&queue);
  return status;
}

crpx_status
crpx_prim (const crpx_graph *graph, crpx_spanning_forest **out)
{ // Grow one minimum spanning tree per connected component, retaining the best crossing edge in a heap.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_status status = graph_check (graph, false);
  if (status != CRPX_OK) return status;
  if (graph->directed) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t size = graph->vertices;
  crpx_spanning_forest *forest = NULL;
  status = forest_new (size, &forest);
  if (status != CRPX_OK) return status;
  crpx_priority_queue queue;
  status = crpx_queue_init (&queue, size);
  if (status != CRPX_OK) { crpx_spanning_forest_free (forest); return status; }
  bool *settled = size ? calloc (size, sizeof (*settled)) : NULL;
  size_t *previous = size ? malloc (size * sizeof (*previous)) : NULL;
  if (size && (!settled || !previous)) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t root = 0; root < size; root++) {
    if (settled[root]) continue;
    forest->components++; previous[root] = SIZE_MAX; crpx_queue_offer (&queue, root, 0);
    while (queue.size) {
      size_t current = crpx_queue_pop (&queue);
      settled[current] = true;
      if (previous[current] != SIZE_MAX)
        forest->edges[forest->edge_count++] = (crpx_graph_edge){previous[current], current, (double)queue.keys[current]};
      for (size_t i = graph->offsets[current]; i < graph->offsets[current + 1]; i++) {
        size_t next = graph->targets[i];
        if (settled[next]) continue;
        if (queue.positions[next] == SIZE_MAX || graph->weights[i] < queue.keys[next]) {
          previous[next] = current; crpx_queue_offer (&queue, next, graph->weights[i]);
        }
      }
    }
  }
  status = forest_total (forest);
done:
  free (settled); free (previous); crpx_queue_free (&queue);
  if (status != CRPX_OK) crpx_spanning_forest_free (forest);
  else *out = forest;
  return status;
}

crpx_status
crpx_kruskal (const crpx_graph *graph, crpx_spanning_forest **out)
{ // Add the cheapest edge joining distinct components; loops and redundant parallel edges are skipped.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_status status = graph_check (graph, false);
  if (status != CRPX_OK) return status;
  if (graph->directed) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t size = graph->vertices;
  if (graph->arcs > SIZE_MAX / sizeof (crpx_graph_edge)) return CRPX_ERROR_OVERFLOW;
  crpx_graph_edge *edges = graph->arcs ? malloc (graph->arcs * sizeof (*edges)) : NULL;
  size_t *parent = size ? malloc (size * sizeof (*parent)) : NULL, *sizes = size ? malloc (size * sizeof (*sizes)) : NULL;
  crpx_spanning_forest *forest = NULL;
  status = forest_new (size, &forest);
  if (status != CRPX_OK) goto done;
  if ((graph->arcs && !edges) || (size && (!parent || !sizes))) { status = CRPX_ERROR_ALLOCATION; goto done; }
  size_t count = 0;
  for (size_t i = 0; i < size; i++) {
    parent[i] = i; sizes[i] = 1;
    for (size_t j = graph->offsets[i]; j < graph->offsets[i + 1]; j++) if (i < graph->targets[j])
      edges[count++] = (crpx_graph_edge){i, graph->targets[j], graph->weights[j]};
  }
  if (count) qsort (edges, count, sizeof (*edges), edge_compare);
  forest->components = size;
  for (size_t i = 0; i < count; i++) {
    size_t first = crpx_component_root (parent, edges[i].first), second = crpx_component_root (parent, edges[i].second);
    if (first == second) continue;
    crpx_component_join (parent, sizes, first, second);
    forest->edges[forest->edge_count++] = edges[i]; forest->components--;
  }
  status = forest_total (forest);
done:
  free (edges); free (parent); free (sizes);
  if (status != CRPX_OK) crpx_spanning_forest_free (forest);
  else *out = forest;
  return status;
}

crpx_status
crpx_prim_dense (size_t vertices, crpx_graph_weight_callback weight, void *data, crpx_spanning_forest **out)
{ // Dense Prim avoids both a matrix and a heap: each newly visited vertex contributes one row of candidate weights.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (vertices && !weight) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_spanning_forest *forest = NULL;
  crpx_status status = forest_new (vertices, &forest);
  if (status != CRPX_OK) return status;
  double *best = vertices ? malloc (vertices * sizeof (*best)) : NULL;
  size_t *previous = vertices ? malloc (vertices * sizeof (*previous)) : NULL;
  bool *settled = vertices ? calloc (vertices, sizeof (*settled)) : NULL;
  if (vertices && (!best || !previous || !settled)) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i = 0; i < vertices; i++) { best[i] = INFINITY; previous[i] = SIZE_MAX; }
  for (size_t position = 0; position < vertices; position++) {
    size_t current = SIZE_MAX;
    for (size_t i = 0; i < vertices; i++) if (!settled[i] && (current == SIZE_MAX || best[i] < best[current])) current = i;
    settled[current] = true;
    if (previous[current] == SIZE_MAX) forest->components++;
    else forest->edges[forest->edge_count++] = (crpx_graph_edge){previous[current], current, best[current]};
    for (size_t i = 0; i < vertices; i++) if (!settled[i]) {
      double candidate;
      status = weight (data, current, i, &candidate);
      if (status != CRPX_OK) goto done;
      if (isnan (candidate) || candidate == -INFINITY) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      if (candidate < best[i]) { best[i] = candidate; previous[i] = current; }
    }
  }
  status = forest_total (forest);
done:
  free (best); free (previous); free (settled);
  if (status != CRPX_OK) crpx_spanning_forest_free (forest);
  else *out = forest;
  return status;
}

static crpx_status
graph_check (const crpx_graph *graph, bool nonnegative)
{
  if (!graph || !graph->offsets || (graph->arcs && (!graph->targets || !graph->weights))) return CRPX_ERROR_INVALID_ARGUMENT;
  if (graph->vertices == SIZE_MAX || graph->vertices > SIZE_MAX / sizeof (long double) ||
      graph->vertices > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  if (graph->offsets[0] || graph->offsets[graph->vertices] != graph->arcs) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < graph->vertices; i++)
    if (graph->offsets[i] > graph->offsets[i + 1] || graph->offsets[i + 1] > graph->arcs) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < graph->arcs; i++)
    if (graph->targets[i] >= graph->vertices || !isfinite (graph->weights[i]) || (nonnegative && graph->weights[i] < 0))
      return CRPX_ERROR_INVALID_ARGUMENT;
  return CRPX_OK;
}

static crpx_status
forest_new (size_t vertices, crpx_spanning_forest **out)
{
  *out = NULL;
  if (vertices > SIZE_MAX / sizeof (size_t) || vertices > SIZE_MAX / sizeof (double) ||
      (vertices && vertices - 1 > SIZE_MAX / sizeof (crpx_graph_edge))) return CRPX_ERROR_OVERFLOW;
  crpx_spanning_forest *forest = calloc (1, sizeof (*forest));
  if (!forest) return CRPX_ERROR_ALLOCATION;
  forest->vertices = vertices;
  forest->edges = vertices > 1 ? malloc ((vertices - 1) * sizeof (*forest->edges)) : NULL;
  if (vertices > 1 && !forest->edges) { free (forest); return CRPX_ERROR_ALLOCATION; }
  *out = forest;
  return CRPX_OK;
}

static crpx_status
forest_total (crpx_spanning_forest *forest)
{
  long double total = 0;
  for (size_t i = 0; i < forest->edge_count; i++) total += forest->edges[i].weight;
  if (!isfinite (total) || fabsl (total) > DBL_MAX) return CRPX_ERROR_OVERFLOW;
  forest->total = (double)total;
  return CRPX_OK;
}

static int
edge_compare (const void *first, const void *second)
{
  const crpx_graph_edge *a = first, *b = second;
  if (a->weight != b->weight) return a->weight < b->weight ? -1 : 1;
  if (a->first != b->first) return a->first < b->first ? -1 : 1;
  return (a->second > b->second) - (a->second < b->second);
}
