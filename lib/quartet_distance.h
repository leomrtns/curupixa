/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file quartet_distance.h
    \brief Compare induced quartet topologies using explicitly corresponding tips, including unresolved quartets.
 */
#ifndef CURUPIXA_QUARTET_DISTANCE_H
#define CURUPIXA_QUARTET_DISTANCE_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Immutable topological path lengths between explicitly ordered tips. The same taxon must occupy the same position in
   indices compared with one another. Each tree may contain additional tips; selected tips must be distinct leaves.
   Branch lengths and the location of the storage root do not affect the score. Public storage is owned by the index;
   treat it as read-only after construction. */
typedef struct {
  size_t taxa;
  uint64_t *distances;
} crpx_quartet_index;
typedef struct {
  uint64_t total, resolved_agreement, resolved_conflict;
  uint64_t first_resolved_only, second_resolved_only, both_unresolved;
} crpx_quartet_counts;
/* Preparation O(nodes log nodes + taxa^2 log nodes), storage O(taxa^2). Multifurcations and unary vertices are
   supported. No inputs are retained. */
crpx_status crpx_quartet_index_new (const crpx_tree *tree, const crpx_node_id *tips, size_t taxa,
                                    crpx_quartet_index **out);
void crpx_quartet_index_free (crpx_quartet_index *index);
/* Four distinct positions: 0 unresolved, 1 ab|cd, 2 ac|bd, 3 ad|bc. */
crpx_status crpx_quartet_topology (const crpx_quartet_index *index, size_t a, size_t b, size_t c, size_t d,
                                   unsigned *topology);
/* Exact enumeration of the taxa-choose-four induced topologies: O(taxa^4), allocation-free. Both unresolved counts as
   agreement; a resolved/unresolved pair counts as disagreement. Counts remain unchanged on error. This reference
   implementation is intended for modest taxon sets, not genome-scale searches. */
crpx_status crpx_quartet_compare (const crpx_quartet_index *first, const crpx_quartet_index *second,
                                  crpx_quartet_counts *counts);
uint64_t crpx_quartet_disagreements (const crpx_quartet_counts *counts);
crpx_status crpx_tree_quartet_distance (const crpx_tree *first, const crpx_node_id *first_tips, const crpx_tree *second,
                                        const crpx_node_id *second_tips, size_t taxa, crpx_quartet_counts *counts);
#ifdef __cplusplus
}
#endif
#endif
