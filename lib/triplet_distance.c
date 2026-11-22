/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "triplet_distance.h"
#include "ancestry_index.h"
#include "topology_counts_internal.h"
#include <stdlib.h>
void
crpx_triplet_index_free (crpx_triplet_index *index)
{
  if (index) {
    free (index->mrca_depth);
    free (index);
  }
}
crpx_status
crpx_triplet_index_new (const crpx_tree *tree, const crpx_node_id *tips, size_t taxa, crpx_triplet_index **out)
{
  /* Store pairwise MRCA depths once; rooted triplet comparisons then require no tree traversal. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || (taxa && !tips)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa && taxa > SIZE_MAX / taxa / sizeof (uint32_t)) return CRPX_ERROR_OVERFLOW;
  crpx_ancestry_index *ancestry = NULL;
  crpx_status status = crpx_ancestry_index_new (tree, &ancestry);
  if (status != CRPX_OK) return status;
  crpx_triplet_index *index = calloc (1, sizeof (*index));
  if (!index) {
    crpx_ancestry_index_free (ancestry);
    return CRPX_ERROR_ALLOCATION;
  }
  index->taxa = taxa;
  index->mrca_depth = taxa ? calloc (taxa * taxa, sizeof (uint32_t)) : NULL;
  if (taxa && !index->mrca_depth) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t i = 0; i < taxa; i++) {
    if (!crpx_ancestry_index_is_leaf (ancestry, tips[i])) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    for (size_t j = 0; j < i; j++) {
      if (tips[i] == tips[j]) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      crpx_node_id ancestor;
      status = crpx_ancestry_index_lca (ancestry, tips[i], tips[j], &ancestor);
      if (status != CRPX_OK) goto done;
      index->mrca_depth[i * taxa + j] = index->mrca_depth[j * taxa + i] =
        crpx_ancestry_index_depth (ancestry, ancestor);
    }
  }
  *out = index;
  index = NULL;
done:
  crpx_triplet_index_free (index);
  crpx_ancestry_index_free (ancestry);
  return status;
}
crpx_status
crpx_triplet_topology (const crpx_triplet_index *index, size_t a, size_t b, size_t c, unsigned *topology)
{
  if (!index || !topology || a >= index->taxa || b >= index->taxa || c >= index->taxa || a == b || a == c || b == c ||
      !index->mrca_depth)
    return CRPX_ERROR_INVALID_ARGUMENT;
  uint32_t ab = index->mrca_depth[a * index->taxa + b], ac = index->mrca_depth[a * index->taxa + c];
  uint32_t bc = index->mrca_depth[b * index->taxa + c];
  *topology = ab > ac && ab > bc ? 1 : ac > ab && ac > bc ? 2 : bc > ab && bc > ac ? 3 : 0;
  return CRPX_OK;
}
crpx_status
crpx_triplet_compare (const crpx_triplet_index *first, const crpx_triplet_index *second, crpx_triplet_counts *counts)
{
  if (!first || !second || !counts || first->taxa != second->taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_triplet_counts result = {0};
  for (size_t a = 0; a < first->taxa; a++)
    for (size_t b = a + 1; b < first->taxa; b++)
      for (size_t c = b + 1; c < first->taxa; c++) {
        unsigned left, right;
        crpx_status status = crpx_triplet_topology (first, a, b, c, &left);
        if (status == CRPX_OK) status = crpx_triplet_topology (second, a, b, c, &right);
        if (status == CRPX_OK) status = crpx_private_topology_counts_add (&result, left, right);
        if (status != CRPX_OK) return status;
      }
  *counts = result;
  return CRPX_OK;
}
crpx_status
crpx_tree_triplet_distance (const crpx_tree *first, const crpx_node_id *first_tips, const crpx_tree *second,
                            const crpx_node_id *second_tips, size_t taxa, crpx_triplet_counts *counts)
{
  if (!counts) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_triplet_index *left = NULL, *right = NULL;
  crpx_status status = crpx_triplet_index_new (first, first_tips, taxa, &left);
  if (status == CRPX_OK) status = crpx_triplet_index_new (second, second_tips, taxa, &right);
  if (status == CRPX_OK) status = crpx_triplet_compare (left, right, counts);
  crpx_triplet_index_free (left);
  crpx_triplet_index_free (right);
  return status;
}
