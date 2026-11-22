/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file paralog_quartets.h
    \brief Prepare and score ASTRAL-Pro quartet classes from rooted, tagged multi-copy gene trees.

    Zhang et al. (2020), ASTRAL-Pro, doi:10.1093/molbev/msaa139, Definitions 1–5. Rooting and event tags are
    supplied by the application.
 */
#ifndef CURUPIXA_PARALOG_QUARTETS_H
#define CURUPIXA_PARALOG_QUARTETS_H
#include "quartet_distance.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t species[4];
  crpx_node_id anchor;
  unsigned topology;
} crpx_paralog_quartet_class;
typedef struct {
  size_t species_taxa, count;
  crpx_paralog_quartet_class *classes;
} crpx_paralog_quartets;
typedef struct {
  uint64_t total, agreement, conflict, unresolved_species;
} crpx_paralog_quartet_score;
/* Exact ASTRAL-Pro per-locus objective for a caller-rooted, caller-tagged binary gene tree (Zhang et al. 2020,
   doi:10.1093/molbev/msaa139, Definitions 1-5). species_of_node maps gene leaves to [0,species_taxa); internal entries
   ignored. duplication_nodes marks internal nodes 0=speciation or 1=duplication. Leaves are ignored. Speciation tags
   must separate disjoint sets of species. Four distinct species and speciation at every triple LCA define eligible
   quartets. Equal species sets and anchor LCA define one class, counted once. Caller tags/root are explicit; no event
   inference or species search occurs. Preparation enumerates O(gene_tips^4) quartets, with logarithmic LCA queries and
   sorting of eligible records. Storage is O(eligible gene-copy quartets) during preparation, then O(classes). Intended
   as an exact reference objective. All input data can be freed after preparation. Public arrays are owned and should
   be treated as read-only. Branch lengths and tip labels are ignored. */
crpx_status crpx_paralog_quartets_new (const crpx_tree *gene, const size_t *species_of_node,
                                      size_t species_taxa, const uint8_t *duplication_nodes,
                                      crpx_paralog_quartets **out);
void crpx_paralog_quartets_free (crpx_paralog_quartets *quartets);
/* Candidate index uses positions matching species_of_node. Allocation-free, O(classes). Unresolved candidate quartets
   contribute no agreement. */
crpx_status crpx_paralog_quartets_score (const crpx_paralog_quartets *quartets,
                                        const crpx_quartet_index *species,
                                        crpx_paralog_quartet_score *score);
#ifdef __cplusplus
}
#endif
#endif
