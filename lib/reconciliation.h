/* Reimplementation of biomcmc-lib reconciliation (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_RECONCILIATION_H
#define CURUPIXA_RECONCILIATION_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_reconciliation crpx_reconciliation;
/* Rooted, strictly bifurcating LCA reconciliation. gene_to_species has one entry
 * per gene node; only gene leaf entries are used and must identify species leaves.
 * Multiple gene copies may map to the same species. The result owns arrays of node
 * IDs/counts and retains no tree pointers. Inputs remain unchanged.
 * Losses count missing descendants below the mapped gene root (the represented
 * family's MRCA); species outside this subtree contribute no losses. This does
 * not optimise over possible gene roots. Extra lineages sum max(crossings-1,0)
 * on species edges using the LCA mapping, excluding the species root stem. */
crpx_status crpx_reconcile_lca(const crpx_tree *gene, const crpx_tree *species,
                              const crpx_node_id *gene_to_species, crpx_reconciliation **out);
void crpx_reconciliation_free(crpx_reconciliation *result);
uint64_t crpx_reconciliation_duplications(const crpx_reconciliation *result);
uint64_t crpx_reconciliation_losses(const crpx_reconciliation *result);
uint64_t crpx_reconciliation_extra_lineages(const crpx_reconciliation *result);
size_t crpx_reconciliation_gene_nodes(const crpx_reconciliation *result);
const crpx_node_id *crpx_reconciliation_mapping(const crpx_reconciliation *result);
const uint8_t *crpx_reconciliation_duplication_nodes(const crpx_reconciliation *result);
const uint64_t *crpx_reconciliation_node_losses(const crpx_reconciliation *result);
const uint64_t *crpx_reconciliation_species_lineages(const crpx_reconciliation *result, size_t *nodes);

#ifdef __cplusplus
}
#endif
#endif
