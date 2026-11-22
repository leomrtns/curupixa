/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Prepare mapped tree comparisons, then call the shared split, dSPR and tanglegram calculations. */
#include "gene_species_distances.h"
#include <stdlib.h>

crpx_status
crpx_gene_species_tanglegram_new (const crpx_tree *gene, const crpx_tree *species, const crpx_node_id *mapping,
                                  crpx_tanglegram_workspace **result)
{ /* Turn each gene-copy mapping into a tip link; retain unmatched species in the second tree. */
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
      links[used++] = (crpx_tanglegram_link){(crpx_node_id) node, mapping[node]};
  status = crpx_tanglegram_workspace_new (gene, species, links, used, result);
  free (links);
  return status;
}

crpx_status
crpx_gene_species_splitsets (const crpx_tree *gene, const crpx_tree *species, const crpx_node_id *mapping,
                             crpx_splitset **gene_splits, crpx_splitset **species_splits)
{ /* Expand the supplied copy map into common split coordinates; the distance algorithms are shared. */
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
  if (!gene_tips || !species_tips) {
    free (gene_tips);
    free (species_tips);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t node = 0; node < crpx_tree_node_count (gene); node++)
    if (crpx_tree_first_child (gene, (crpx_node_id) node) == CRPX_NODE_NONE) {
      gene_tips[count] = (crpx_node_id) node;
      species_tips[count++] = mapping[node];
    }
  status = crpx_splitset_from_tree (gene, gene_tips, taxa, &first);
  if (status == CRPX_OK) status = crpx_splitset_from_tree (species, species_tips, taxa, &second);
  free (gene_tips);
  free (species_tips);
  if (status != CRPX_OK) {
    crpx_splitset_free (first);
    crpx_splitset_free (second);
    return status;
  }
  *gene_splits = first;
  *species_splits = second;
  return CRPX_OK;
}

crpx_status
crpx_gene_species_split_distances (const crpx_tree *gene, const crpx_tree *species, const crpx_node_id *mapping,
                                   size_t *rf, uint64_t *assignment_cost)
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
crpx_gene_species_spr_distance (const crpx_tree *gene, const crpx_tree *species, const crpx_node_id *mapping,
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

void
crpx_gene_species_comparison_free (crpx_gene_species_comparison *comparison)
{
  if (!comparison) return;
  crpx_tree_comparison_free (comparison->gene);
  crpx_tree_comparison_free (comparison->species);
  free (comparison);
}
crpx_status
crpx_gene_species_comparison_new (const crpx_tree *gene, const crpx_tree *species, const crpx_node_id *mapping,
                                  bool rooted, bool with_lengths, crpx_gene_species_comparison **out)
{ /* All expansion happens once here. Repeated distance calls use the unchanged numerical kernels. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!gene || !species || !mapping) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t copies = crpx_tree_leaf_count (gene), nodes = crpx_tree_node_count (species);
  if (!copies || !nodes) return CRPX_ERROR_INVALID_ARGUMENT;
  if (copies > SIZE_MAX / sizeof (crpx_node_id) || nodes > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  crpx_gene_species_comparison *comparison = calloc (1, sizeof (*comparison));
  size_t *counts = calloc (nodes, sizeof (size_t));
  crpx_node_id *first = malloc (copies * sizeof (*first)), *second = malloc (copies * sizeof (*second));
  crpx_tree *expanded = NULL;
  crpx_tree_builder *builder = NULL;
  crpx_status status = CRPX_OK;
  if (!comparison || !counts || !first || !second) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  size_t used = 0;
  for (size_t node = 0; node < crpx_tree_node_count (gene); node++)
    if (crpx_tree_first_child (gene, (crpx_node_id) node) == CRPX_NODE_NONE) {
      crpx_node_id target = mapping[node];
      if (target >= nodes || crpx_tree_first_child (species, target) != CRPX_NODE_NONE) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      if (!counts[target]++) comparison->represented_species++;
      first[used] = (crpx_node_id) node;
      second[used++] = target;
    }
  comparison->gene_copies = copies;
  comparison->total_species = crpx_tree_leaf_count (species);
  status = crpx_tree_clone (species, &expanded);
  if (status == CRPX_OK) status = crpx_tree_builder_new (&builder);
  if (status != CRPX_OK) goto done;
  for (size_t i = 0; i < copies; i++)
    if (counts[second[i]] > 1) {
      crpx_node_id copy;
      status = crpx_tree_builder_append (builder, expanded, second[i], &copy);
      if (status == CRPX_OK) status = crpx_tree_set_branch_length (expanded, copy, 0);
      if (status != CRPX_OK) goto done;
      second[i] = copy;
    }
  status = crpx_tree_comparison_new (gene, first, copies, rooted, with_lengths, &comparison->gene);
  if (status == CRPX_OK)
    status = crpx_tree_comparison_new (expanded, second, copies, rooted, with_lengths, &comparison->species);
  if (status == CRPX_OK) {
    *out = comparison;
    comparison = NULL;
  }
done:
  crpx_gene_species_comparison_free (comparison);
  crpx_tree_builder_free (builder);
  crpx_tree_free (expanded);
  free (counts);
  free (first);
  free (second);
  return status;
}
