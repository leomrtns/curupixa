/* Derived from guenomu/genefam-dist topology_mrca.
 * Copyright (C) 2009-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_RECONCILIATION_ROOTS_H
#define CURUPIXA_RECONCILIATION_ROOTS_H
#include "ancestry_index.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { CRPX_RECONCILE_ALL_SPECIES, CRPX_RECONCILE_REPRESENTED_SPECIES } crpx_reconciliation_scope;
typedef struct {
  crpx_node_id first,second; /* endpoints of an unrooted gene edge, original IDs */
  uint64_t duplications,losses,extra_lineages;
} crpx_root_reconciliation_score;
typedef struct crpx_root_reconciliations crpx_root_reconciliations;
/* Score every edge-rooting of a strictly binary gene tree against a binary species
 * snapshot, using two directional passes. The artificial binary gene root is
 * suppressed; its two children identify the single bridged edge. Branch lengths
 * and labels are ignored, and neither tree nor index is modified or retained.
 * Explicit map has gene_nodes entries; only leaf entries are read.
 * ALL_SPECIES counts losses on the full species tree, within the family's MRCA.
 * REPRESENTED_SPECIES suppresses unrepresented side branches (legacy convention).
 * Extra lineages always use the represented species tree and the legacy identity
 * losses_represented-2*duplications+2*(gene_leaves-distinct_species).
 * Each objective has its own minimum: do not combine minima into a fictitious
 * single rooting. Array order is deterministic by endpoint IDs. A one-tip gene
 * has one degenerate candidate with both endpoints NONE and zero scores.
 * O(species_nodes + gene_nodes*log(species_nodes)) work; O(species_nodes+gene_nodes)
 * extra memory, in addition to the reusable ancestry index. */
crpx_status crpx_reconcile_all_roots(const crpx_tree *gene,const crpx_ancestry_index *species,
                                    const crpx_node_id *gene_to_species,crpx_reconciliation_scope scope,
                                    crpx_root_reconciliations **out);
void crpx_root_reconciliations_free(crpx_root_reconciliations *result);
const crpx_root_reconciliation_score *crpx_root_reconciliations_scores(const crpx_root_reconciliations *result,
                                                                     size_t *count);
/* Indices into scores(), with first encountered tie. */
crpx_status crpx_root_reconciliations_minima(const crpx_root_reconciliations *result,
                                           size_t *duplications,size_t *losses,size_t *extra_lineages);
#ifdef __cplusplus
}
#endif
#endif
