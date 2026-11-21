/* SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file gene_species_distances.h
    \brief Prepare mapped tree comparisons, then call the shared split, dSPR and tanglegram calculations.
 */
#ifndef CURUPIXA_GENE_SPECIES_DISTANCES_H
#define CURUPIXA_GENE_SPECIES_DISTANCES_H
#include "spr_distance.h"
#include "tanglegram.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Prepare comparable splitsets once. Bit order follows gene leaf node IDs. Mapping is indexed by gene node; each gene
   leaf identifies a species leaf. Several copies may identify the same species. Unrepresented species add no bits. The
   results own their data and retain no input pointers. Use the ordinary splitset comparison/SPR functions repeatedly
   without remapping either tree. Rebuild the affected splitset after a topology or mapping change. */
crpx_status crpx_gene_species_splitsets (const crpx_tree *gene, const crpx_tree *species,
                                        const crpx_node_id *mapping,
                                        crpx_splitset **gene_splits, crpx_splitset **species_splits);
crpx_status crpx_gene_species_split_distances (const crpx_tree *gene, const crpx_tree *species,
                                             const crpx_node_id *mapping,
                                             size_t *rf, uint64_t *assignment_cost);
crpx_status crpx_gene_species_spr_distance (const crpx_tree *gene, const crpx_tree *species,
                                          const crpx_node_id *mapping,
                                          crpx_spr_workspace *workspace, crpx_spr_result *result);
/* Prepare one link per gene tip to its mapped species tip. Species without copies remain unlinked. Calls the general
   unequal-size tanglegram kernel. */
crpx_status crpx_gene_species_tanglegram_new (const crpx_tree *gene, const crpx_tree *species,
                                            const crpx_node_id *mapping,
                                            crpx_tanglegram_workspace **result);
#ifdef __cplusplus
}
#endif
#endif
