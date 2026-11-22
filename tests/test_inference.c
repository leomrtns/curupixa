/* SPDX-License-Identifier: MIT */
#include "inference.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static crpx_node_id
leaf (const crpx_tree *tree, const char *name)
{
  for (size_t i = 0; i < crpx_tree_node_count (tree); i++) {
    const char *label = crpx_tree_label (tree, (crpx_node_id) i);
    if (label && !strcmp (label, name)) return (crpx_node_id) i;
  }
  return CRPX_NODE_NONE;
}

int main (void)
{
  bool accepted = false;
  OK (crpx_metropolis_accept (1, 0, log (0.9), &accepted)); assert (accepted);
  OK (crpx_metropolis_accept (-2, 0, log (0.1), &accepted)); assert (accepted);
  OK (crpx_metropolis_accept (-2, 0, log (0.9), &accepted)); assert (!accepted);
  OK (crpx_metropolis_accept (-2, 3, log (0.9), &accepted)); assert (accepted);
  OK (crpx_metropolis_accept (-INFINITY, 0, -INFINITY, &accepted)); assert (!accepted);
  OK (crpx_metropolis_accept (INFINITY, 0, log (0.9), &accepted)); assert (accepted);
  assert (crpx_metropolis_accept (INFINITY, -INFINITY, log (0.5), &accepted) == CRPX_ERROR_INVALID_ARGUMENT);
  assert (accepted);
  assert (crpx_metropolis_accept (0, 0, 0, &accepted) == CRPX_ERROR_INVALID_ARGUMENT);
  assert (crpx_metropolis_accept (NAN, 0, -1, &accepted) == CRPX_ERROR_INVALID_ARGUMENT);

  crpx_tree *gene, *species;
  OK (crpx_tree_from_newick ("(A,B);", &gene, NULL));
  OK (crpx_tree_from_newick ("(A,B);", &species, NULL));
  size_t nodes = crpx_tree_node_count (gene);
  crpx_node_id *map = malloc (nodes * sizeof (*map)); assert (map);
  for (size_t i = 0; i < nodes; i++) map[i] = CRPX_NODE_NONE;
  map[leaf (gene, "A")] = leaf (species, "A");
  map[leaf (gene, "B")] = leaf (species, "B");
  crpx_dtl_species *prepared;
  crpx_dtl_workspace *workspace;
  crpx_dtl_costs costs = {2, 3, 1};
  crpx_dtl_history_probability result;
  OK (crpx_dtl_species_new (species, NULL, 0, &prepared));
  OK (crpx_dtl_workspace_new (&workspace));
  OK (crpx_dtl_minimum_log_probability (gene, prepared, map, &costs, workspace, &result));
  assert (result.minimum_cost == 0 && result.log_partition >= 0 && result.log_probability <= 0);
  assert (fabs (result.log_probability + result.minimum_cost + result.log_partition) < 1e-12);
  crpx_dtl_mapping *mapping = malloc (nodes * sizeof (*mapping)); assert (mapping);
  crpx_dtl_counts counts;
  OK (crpx_dtl_traceback (workspace, mapping, nodes, &counts));
  assert (!counts.duplications && !counts.losses && !counts.transfers);
  double partition;
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_LOG_PARTITION, workspace, &partition));
  assert (fabs (partition - result.log_partition) < 1e-12);
  costs.transfer = 0;
  OK (crpx_dtl_minimum_log_probability (gene, prepared, map, &costs, workspace, &result));
  assert (result.log_partition >= partition);
  free (map); free (mapping);
  crpx_dtl_workspace_free (workspace); crpx_dtl_species_free (prepared);
  crpx_tree_free (gene); crpx_tree_free (species);
  puts ("inference tests passed");
  return 0;
}
