/* Reimplementation of biomcmc-lib splitset distances (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_SPLITSET_DISTANCES_H
#define CURUPIXA_SPLITSET_DISTANCES_H
#include "tree.h"
#include "bipartition.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_splitset crpx_splitset;
/* Bit i represents tip_nodes[i]. Repeated node IDs expand a species leaf into
 * multiple gene copies. Unrepresented leaves contribute no bits. Extract unique
 * unrooted nontrivial splits, including expanded multi-copy terminal groups.
 * Taxon bit order must match between splitsets being compared. */
crpx_status crpx_splitset_from_tree(const crpx_tree *tree, const crpx_node_id *tip_nodes, size_t taxa,
                                   crpx_splitset **out);
void crpx_splitset_free(crpx_splitset *splits);
size_t crpx_splitset_count(const crpx_splitset *splits);
const struct crpx_bipartition_struct *crpx_splitset_get(const crpx_splitset *splits, size_t index);
/* RF is the symmetric difference. Assignment cost matches disagreeing splits by
 * min(|a XOR b|, taxa-|a XOR b|). An unmatched split costs its smaller-side size.
 * This positive dummy-edge penalty intentionally replaces biomcmc's zero penalty.
 * Either output may be NULL. RF-only uses expected linear hash matching in the
 * number of split words; assignment additionally builds the quadratic cost matrix. */
crpx_status crpx_splitset_compare(const crpx_splitset *a, const crpx_splitset *b,
                                 size_t *rf, uint64_t *assignment_cost);
/* gene_to_species is indexed by gene node ID; only gene leaves are used. */
crpx_status crpx_gene_species_split_distances(const crpx_tree *gene, const crpx_tree *species,
                                             const crpx_node_id *gene_to_species,
                                             size_t *rf, uint64_t *assignment_cost);
#ifdef __cplusplus
}
#endif
#endif
