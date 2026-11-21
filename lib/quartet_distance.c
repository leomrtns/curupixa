/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "quartet_distance.h"
#include "ancestry_index.h"
#include <stdlib.h>

void
crpx_quartet_index_free (crpx_quartet_index *index)
{
  if (!index) return;
  free (index->distances);
  free (index);
}

crpx_status
crpx_quartet_index_new (const crpx_tree *tree, const crpx_node_id *tips,
                        size_t taxa, crpx_quartet_index **out)
{
  crpx_ancestry_index *ancestry = NULL;
  crpx_quartet_index *index;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || (taxa && !tips)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa && taxa > SIZE_MAX / taxa / sizeof (uint64_t)) return CRPX_ERROR_OVERFLOW;
  status = crpx_ancestry_index_new (tree, &ancestry);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < taxa; i++) {
    if (!crpx_ancestry_index_is_leaf (ancestry, tips[i])) {
      crpx_ancestry_index_free (ancestry);
      return CRPX_ERROR_INVALID_ARGUMENT;
    }
    for (size_t j = 0; j < i; j++) if (tips[i] == tips[j]) {
      crpx_ancestry_index_free (ancestry);
      return CRPX_ERROR_INVALID_ARGUMENT;
    }
  }
  index = calloc (1, sizeof (*index));
  if (!index) { crpx_ancestry_index_free (ancestry); return CRPX_ERROR_ALLOCATION; }
  index->taxa = taxa;
  index->distances = taxa ? calloc (taxa * taxa, sizeof (*index->distances)) : NULL;
  if (taxa && !index->distances) {
    crpx_ancestry_index_free (ancestry);
    crpx_quartet_index_free (index);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < taxa; i++) for (size_t j = 0; j < i; j++) {
    crpx_node_id ancestor;
    status = crpx_ancestry_index_lca (ancestry, tips[i], tips[j], &ancestor);
    if (status != CRPX_OK) {
      crpx_ancestry_index_free (ancestry);
      crpx_quartet_index_free (index);
      return status;
    }
    uint64_t distance = (uint64_t) crpx_ancestry_index_depth (ancestry, tips[i]) +
      crpx_ancestry_index_depth (ancestry, tips[j]) - 2 * (uint64_t) crpx_ancestry_index_depth (ancestry, ancestor);
    index->distances[i * taxa + j] = index->distances[j * taxa + i] = distance;
  }
  crpx_ancestry_index_free (ancestry);
  *out = index;
  return CRPX_OK;
}

static unsigned
quartet_topology (const crpx_quartet_index *index, size_t a, size_t b, size_t c, size_t d)
{
  size_t n = index->taxa;
  const uint64_t *distance = index->distances;
  uint64_t ab_cd = distance[a * n + b] + distance[c * n + d];
  uint64_t ac_bd = distance[a * n + c] + distance[b * n + d];
  uint64_t ad_bc = distance[a * n + d] + distance[b * n + c];
  if (ab_cd < ac_bd && ab_cd < ad_bc) return 1;
  if (ac_bd < ab_cd && ac_bd < ad_bc) return 2;
  if (ad_bc < ab_cd && ad_bc < ac_bd) return 3;
  return 0;
}

crpx_status
crpx_quartet_topology (const crpx_quartet_index *index, size_t a, size_t b,
                       size_t c, size_t d, unsigned *topology)
{
  if (!index || !topology || !index->distances || a >= index->taxa || b >= index->taxa ||
      c >= index->taxa || d >= index->taxa || a == b || a == c || a == d || b == c || b == d || c == d)
    return CRPX_ERROR_INVALID_ARGUMENT;
  *topology = quartet_topology (index, a, b, c, d);
  return CRPX_OK;
}

crpx_status
crpx_quartet_compare (const crpx_quartet_index *first, const crpx_quartet_index *second,
                      crpx_quartet_counts *counts)
{
  crpx_quartet_counts result = {0};
  if (!first || !second || !counts || first->taxa != second->taxa ||
      (first->taxa && (!first->distances || !second->distances))) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = first->taxa;
  for (size_t a = 0; a < n; a++) for (size_t b = a + 1; b < n; b++)
    for (size_t c = b + 1; c < n; c++) for (size_t d = c + 1; d < n; d++) {
      if (result.total == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
      result.total++;
      unsigned left = quartet_topology (first, a, b, c, d), right = quartet_topology (second, a, b, c, d);
      if (!left && !right) result.both_unresolved++;
      else if (!left) result.second_resolved_only++;
      else if (!right) result.first_resolved_only++;
      else if (left == right) result.resolved_agreement++;
      else result.resolved_conflict++;
    }
  *counts = result;
  return CRPX_OK;
}

uint64_t
crpx_quartet_disagreements (const crpx_quartet_counts *counts)
{
  return counts ? counts->resolved_conflict + counts->first_resolved_only + counts->second_resolved_only : 0;
}

crpx_status
crpx_tree_quartet_distance (const crpx_tree *first, const crpx_node_id *first_tips,
                            const crpx_tree *second, const crpx_node_id *second_tips,
                            size_t taxa, crpx_quartet_counts *counts)
{
  crpx_quartet_index *left = NULL, *right = NULL;
  if (!counts) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_quartet_index_new (first, first_tips, taxa, &left);
  if (status == CRPX_OK) status = crpx_quartet_index_new (second, second_tips, taxa, &right);
  if (status == CRPX_OK) status = crpx_quartet_compare (left, right, counts);
  crpx_quartet_index_free (left);
  crpx_quartet_index_free (right);
  return status;
}
