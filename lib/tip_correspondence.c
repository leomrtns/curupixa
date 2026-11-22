/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "tip_correspondence.h"
#include "hashtable.h"
#include <stdlib.h>
void
crpx_tip_correspondence_free (crpx_tip_correspondence *correspondence)
{
  if (correspondence) {
    free (correspondence->first);
    free (correspondence->second);
    free (correspondence);
  }
}
crpx_status
crpx_tip_correspondence_new (const crpx_tree *first, const crpx_tree *second, crpx_tip_correspondence **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!first || !second) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa = crpx_tree_leaf_count (first);
  if (!taxa || taxa != crpx_tree_leaf_count (second)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX / sizeof (crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = crpx_tree_validate_policy (first, CRPX_TREE_REQUIRE_UNIQUE_TIPS);
  if (status == CRPX_OK) status = crpx_tree_validate_policy (second, CRPX_TREE_REQUIRE_UNIQUE_TIPS);
  if (status != CRPX_OK) return status;
  crpx_tip_correspondence *correspondence = calloc (1, sizeof (*correspondence));
  crpx_hashtable_t labels = NULL;
  if (!correspondence) return CRPX_ERROR_ALLOCATION;
  correspondence->taxa = taxa;
  correspondence->first = malloc (taxa * sizeof (crpx_node_id));
  correspondence->second = malloc (taxa * sizeof (crpx_node_id));
  if (!correspondence->first || !correspondence->second) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  status = crpx_hashtable_new (taxa, 0, &labels);
  if (status != CRPX_OK) goto done;
  for (size_t node = 0; node < crpx_tree_node_count (second); node++)
    if (crpx_tree_first_child (second, (crpx_node_id) node) == CRPX_NODE_NONE) {
      status = crpx_hashtable_put_string (labels, crpx_tree_label (second, (crpx_node_id) node), node, NULL);
      if (status != CRPX_OK) goto done;
    }
  size_t used = 0;
  for (size_t node = 0; node < crpx_tree_node_count (first); node++)
    if (crpx_tree_first_child (first, (crpx_node_id) node) == CRPX_NODE_NONE) {
      uint64_t target;
      status = crpx_hashtable_get_string (labels, crpx_tree_label (first, (crpx_node_id) node), &target);
      if (status == CRPX_ERROR_NOT_FOUND) status = CRPX_ERROR_INVALID_ARGUMENT;
      if (status != CRPX_OK) goto done;
      correspondence->first[used] = (crpx_node_id) node;
      correspondence->second[used++] = (crpx_node_id) target;
    }
  *out = correspondence;
  correspondence = NULL;
done:
  crpx_hashtable_free (labels);
  crpx_tip_correspondence_free (correspondence);
  return status;
}
