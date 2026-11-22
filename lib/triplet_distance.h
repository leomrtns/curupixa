/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file triplet_distance.h
    \brief Compare rooted triplets with the same count and unresolved conventions as quartets. */
#ifndef CURUPIXA_TRIPLET_DISTANCE_H
#define CURUPIXA_TRIPLET_DISTANCE_H
#include "quartet_distance.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
  size_t taxa;
  uint32_t *mrca_depth;
} crpx_triplet_index;
typedef crpx_quartet_counts crpx_triplet_counts;
/* Explicit distinct tip correspondence; selected subsets are permitted but never
   inferred by intersecting labels. Root orientation matters; lengths do not.
   Preparation O(nodes log nodes + taxa^2 log nodes); comparison O(taxa^3), no allocation. */
crpx_status crpx_triplet_index_new (const crpx_tree *tree, const crpx_node_id *tips, size_t taxa,
                                    crpx_triplet_index **out);
void crpx_triplet_index_free (crpx_triplet_index *index);
/* Topology 0=unresolved, 1=ab|c, 2=ac|b, 3=bc|a. */
crpx_status crpx_triplet_topology (const crpx_triplet_index *index, size_t a, size_t b, size_t c, unsigned *topology);
crpx_status crpx_triplet_compare (const crpx_triplet_index *first, const crpx_triplet_index *second,
                                  crpx_triplet_counts *counts);
crpx_status crpx_tree_triplet_distance (const crpx_tree *first, const crpx_node_id *first_tips, const crpx_tree *second,
                                        const crpx_node_id *second_tips, size_t taxa, crpx_triplet_counts *counts);
#ifdef __cplusplus
}
#endif
#endif
