/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file graph.h
    Shortest paths and minimum spanning forests for arbitrary weighted graphs and supplied pairwise weights.
    No tree, sequence or clustering representation is required. */
#ifndef CURUPIXA_GRAPH_H
#define CURUPIXA_GRAPH_H
#include <stdbool.h>
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { size_t first, second; double weight; } crpx_graph_edge;
/* Construct through graph_new; arrays belong to the graph and remain unchanged during calculations. Undirected
   edges are stored in both directions, except loops. Parallel edges and zero weights are permitted. */
typedef struct {
  size_t vertices, arcs;
  size_t *offsets, *targets;
  double *weights;
  bool directed;
} crpx_graph;

typedef struct {
  size_t vertices, edge_count, components;
  crpx_graph_edge *edges;
  double total;
} crpx_spanning_forest;

/* Copies finite edge weights, including negative values. Vertex IDs are 0..vertices-1. */
crpx_status crpx_graph_new (size_t vertices, const crpx_graph_edge *edges, size_t count, bool directed, crpx_graph **out);
void crpx_graph_free (crpx_graph *graph);
void crpx_spanning_forest_free (crpx_spanning_forest *forest);

/* Dijkstra (1959): non-negative weights, O((vertices+arcs) log vertices) time and O(vertices) scratch storage.
   Outputs use vertex IDs; unreachable distances are INFINITY and missing predecessors SIZE_MAX (also at the source).
   Outputs remain unchanged on error, including unrepresentable finite shortest distances. */
crpx_status crpx_dijkstra (const crpx_graph *graph, size_t source, double *distances, size_t *predecessors, size_t capacity);
/* Prim (1957) and Kruskal (1956): undirected graphs, including negative weights and disconnected components.
   A minimum spanning forest contains vertices-components edges. Prim uses a heap; Kruskal sorts the edges. */
crpx_status crpx_prim (const crpx_graph *graph, crpx_spanning_forest **out);
crpx_status crpx_kruskal (const crpx_graph *graph, crpx_spanning_forest **out);

/* Dense Prim queries only distinct pairs; weights must be symmetric, finite or +INFINITY (no edge).
   It uses O(vertices) scratch storage and O(vertices^2) queries, without constructing a distance matrix. */
typedef crpx_status (*crpx_graph_weight_callback) (void *data, size_t first, size_t second, double *weight);
crpx_status crpx_prim_dense (size_t vertices, crpx_graph_weight_callback weight, void *data, crpx_spanning_forest **out);

#ifdef __cplusplus
}
#endif
#endif
