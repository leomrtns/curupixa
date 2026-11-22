/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file hdbscan.h
    HDBSCAN* from general dissimilarities or a supplied mutual-reachability minimum spanning tree.
    Condensed edges and selected cluster stabilities remain available for downstream analysis. */
#ifndef CURUPIXA_HDBSCAN_H
#define CURUPIXA_HDBSCAN_H
#include <stdbool.h>
#include <stddef.h>
#include "distance_generator.h"
#include "graph.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CRPX_HDBSCAN_EOM = 0, CRPX_HDBSCAN_LEAF = 1 } crpx_hdbscan_selection;
typedef struct { size_t parent, child, size; double lambda; } crpx_hdbscan_edge;
typedef struct {
  size_t samples, clusters, edge_count;
  size_t *labels, *selected_clusters;
  double *probabilities, *stabilities;
  crpx_hdbscan_edge *hierarchy;
} crpx_hdbscan;

/* min_points includes the query itself, as in OPTICS. min_cluster_size lies in 2..samples.
   EOM selects a non-overlapping set maximising stability; LEAF selects leaves of the condensed hierarchy.
   EOM root selection requires allow_single_cluster=true; LEAF selects non-root leaves.
   Noise labels are SIZE_MAX, with membership zero.
   labels index selected_clusters/stabilities. Condensed point IDs are 0..samples-1; cluster IDs start at samples.
   Distances must be finite, symmetric and non-negative. O(n^2) distance queries and O(n) extra storage;
   the optional generator cache has its own storage cost. No coordinate index or approximate MST is assumed.
   Campello et al. (2015), doi:10.1145/2733381. */
crpx_status crpx_hdbscan_run (crpx_distance_generator_t generator, size_t component, size_t min_points,
                             size_t min_cluster_size, crpx_hdbscan_selection selection, bool allow_single_cluster,
                             crpx_hdbscan **out);
/* Supply an already computed, connected mutual-reachability MST to reuse an external neighbour index or cache.
   The function checks tree structure and non-negative weights; it cannot verify that supplied weights are an MST
   of max(core[first],core[second],distance). Input storage is borrowed only during this call. */
crpx_status crpx_hdbscan_from_mst (const crpx_spanning_forest *mst, size_t min_cluster_size,
                                  crpx_hdbscan_selection selection, bool allow_single_cluster, crpx_hdbscan **out);
void crpx_hdbscan_free (crpx_hdbscan *result);
#ifdef __cplusplus
}
#endif
#endif
