/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "gene_species_distances.h"
#include <stdlib.h>

crpx_status
crpx_gene_species_tanglegram_new (const crpx_tree *gene, const crpx_tree *species,
                                 const crpx_node_id *mapping, crpx_tanglegram_workspace **result)
{
  size_t taxa = crpx_tree_leaf_count (gene), used = 0;
  crpx_tanglegram_link *links;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!gene || !species || !mapping || !taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX / sizeof (*links)) return CRPX_ERROR_OVERFLOW;
  links = malloc (taxa * sizeof (*links));
  if (!links) return CRPX_ERROR_ALLOCATION;
  for (size_t node = 0; node < crpx_tree_node_count (gene); node++)
    if (crpx_tree_first_child (gene, (crpx_node_id) node) == CRPX_NODE_NONE)
      links[used++] = (crpx_tanglegram_link) {(crpx_node_id) node, mapping[node]};
  status = crpx_tanglegram_workspace_new (gene, species, links, used, result);
  free (links);
  return status;
}

crpx_status
crpx_gene_species_splitsets (const crpx_tree *gene, const crpx_tree *species,
                            const crpx_node_id *mapping,
                            crpx_splitset **gene_splits, crpx_splitset **species_splits)
{
  size_t taxa = crpx_tree_leaf_count (gene), count = 0;
  crpx_node_id *gene_tips, *species_tips;
  crpx_splitset *first = NULL, *second = NULL;
  crpx_status status;
  if (!gene_splits || !species_splits || gene_splits == species_splits) return CRPX_ERROR_INVALID_ARGUMENT;
  *gene_splits = NULL;
  *species_splits = NULL;
  if (!gene || !species || !mapping || !taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX / sizeof (*gene_tips)) return CRPX_ERROR_OVERFLOW;
  gene_tips = malloc (taxa * sizeof (*gene_tips));
  species_tips = malloc (taxa * sizeof (*species_tips));
  if (!gene_tips || !species_tips) { free (gene_tips); free (species_tips); return CRPX_ERROR_ALLOCATION; }
  for (size_t node = 0; node < crpx_tree_node_count (gene); node++)
    if (crpx_tree_first_child (gene, (crpx_node_id) node) == CRPX_NODE_NONE) {
      gene_tips[count] = (crpx_node_id) node;
      species_tips[count++] = mapping[node];
    }
  status = crpx_splitset_from_tree (gene, gene_tips, taxa, &first);
  if (status == CRPX_OK) status = crpx_splitset_from_tree (species, species_tips, taxa, &second);
  free (gene_tips);
  free (species_tips);
  if (status != CRPX_OK) { crpx_splitset_free (first); crpx_splitset_free (second); return status; }
  *gene_splits = first;
  *species_splits = second;
  return CRPX_OK;
}

crpx_status
crpx_gene_species_split_distances (const crpx_tree *gene, const crpx_tree *species,
                                  const crpx_node_id *mapping, size_t *rf, uint64_t *assignment_cost)
{
  crpx_splitset *first = NULL, *second = NULL;
  crpx_status status;
  if (!rf && !assignment_cost) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_gene_species_splitsets (gene, species, mapping, &first, &second);
  if (status == CRPX_OK) status = crpx_splitset_compare (first, second, rf, assignment_cost);
  crpx_splitset_free (first);
  crpx_splitset_free (second);
  return status;
}

crpx_status
crpx_gene_species_spr_distance (const crpx_tree *gene, const crpx_tree *species,
                               const crpx_node_id *mapping,
                               crpx_spr_workspace *workspace, crpx_spr_result *result)
{
  crpx_splitset *first = NULL, *second = NULL;
  crpx_status status;
  if (!workspace || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_gene_species_splitsets (gene, species, mapping, &first, &second);
  if (status == CRPX_OK) status = crpx_splitset_spr_distance (first, second, workspace, result);
  crpx_splitset_free (first);
  crpx_splitset_free (second);
  return status;
}
