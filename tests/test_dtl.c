/* SPDX-License-Identifier: MIT */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "dtl.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static void
map_tips (const crpx_tree *gene, const crpx_tree *species, crpx_node_id *mapping)
{
  for (size_t node = 0; node < crpx_tree_node_count (gene); node++) {
    mapping[node] = CRPX_NODE_NONE;
    if (crpx_tree_first_child (gene, (crpx_node_id) node) != CRPX_NODE_NONE) continue;
    const char *label = crpx_tree_label (gene, (crpx_node_id) node);
    size_t length = strcspn (label, "|");
    for (size_t tip = 0; tip < crpx_tree_node_count (species); tip++) {
      if (crpx_tree_first_child (species, (crpx_node_id) tip) != CRPX_NODE_NONE) continue;
      const char *species_label = crpx_tree_label (species, (crpx_node_id) tip);
      if (strlen (species_label) == length && !strncmp (label, species_label, length))
        mapping[node] = (crpx_node_id) tip;
    }
    assert (mapping[node] != CRPX_NODE_NONE);
  }
}

static void
check_two_species (bool dated)
{
  crpx_tree *gene, *species;
  OK (crpx_tree_from_newick ("(A|1,B|1);", &gene, NULL));
  OK (crpx_tree_from_newick ("(A,B);", &species, NULL));
  crpx_node_id map[3];
  map_tips (gene, species, map);
  double ages[3];
  for (size_t node = 0; node < 3; node++)
    ages[node] = crpx_tree_first_child (species, (crpx_node_id) node) == CRPX_NODE_NONE ? 0 : 1;
  crpx_dtl_species *prepared;
  OK (crpx_dtl_species_new (species, dated ? ages : NULL, 2, &prepared));
  crpx_dtl_workspace *workspace;
  OK (crpx_dtl_workspace_new (&workspace));
  crpx_dtl_costs costs = {3, 2, 1};
  double value;
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value));
  assert (value == 0);
  size_t storage = crpx_dtl_workspace_storage_size (workspace);
  crpx_dtl_mapping mapping[3];
  crpx_dtl_counts counts;
  OK (crpx_dtl_traceback (workspace, mapping, 3, &counts));
  assert (!counts.duplications && !counts.transfers && !counts.losses);
  assert (mapping[crpx_tree_root (gene)].event == CRPX_DTL_SPECIATION);
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_LOG_PARTITION, workspace, &value));
  /* Enumerate the four observed histories: speciation; loss then A->B transfer;
   * loss then B->A transfer; root-stem duplication with two subsequent losses. */
  double partition = 1 + 2 * exp (-costs.loss - costs.transfer) + exp (-costs.duplication - 2 * costs.loss);
  assert (fabs (value - log (partition)) < 1e-12);
  assert (crpx_dtl_traceback (workspace, mapping, 3, &counts) == CRPX_ERROR_INVALID_ARGUMENT);
  assert (storage == crpx_dtl_workspace_storage_size (workspace));
  costs.transfer = 0;
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value));
  assert (value == 0);
  map[1] = CRPX_NODE_NONE;
  value = 123;
  assert (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value) ==
          CRPX_ERROR_INVALID_ARGUMENT && value == 123);
  crpx_dtl_workspace_free (workspace);
  crpx_dtl_species_free (prepared);
  crpx_tree_free (gene);
  crpx_tree_free (species);
}

static void
check_paralogs_and_losses (void)
{
  crpx_tree *gene, *species;
  OK (crpx_tree_from_newick ("(A|1,A|2);", &gene, NULL));
  OK (crpx_tree_from_newick ("(A,B);", &species, NULL));
  crpx_node_id map[3];
  map_tips (gene, species, map);
  crpx_dtl_species *prepared;
  OK (crpx_dtl_species_new (species, NULL, 0, &prepared));
  crpx_dtl_workspace *workspace;
  OK (crpx_dtl_workspace_new (&workspace));
  crpx_dtl_costs costs = {3, 2, 1};
  double value;
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value));
  assert (value == 4);
  crpx_dtl_mapping mapping[3];
  crpx_dtl_counts counts;
  OK (crpx_dtl_traceback (workspace, mapping, 3, &counts));
  assert (counts.duplications == 1 && counts.losses == 1 && counts.transfers == 0);
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_LOG_PARTITION, workspace, &value));
  assert (fabs (value - log (exp (-4) + exp (-5))) < 1e-12);
  crpx_dtl_species_free (prepared);
  crpx_tree_free (species);
  OK (crpx_tree_from_newick ("A;", &species, NULL));
  map_tips (gene, species, map);
  OK (crpx_dtl_species_new (species, NULL, 0, &prepared));
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value));
  assert (value == 3);
  OK (crpx_dtl_traceback (workspace, mapping, 3, &counts));
  assert (counts.duplications == 1 && !counts.losses && !counts.transfers);
  crpx_dtl_workspace_clear (workspace);
  crpx_dtl_workspace_free (workspace);
  crpx_dtl_species_free (prepared);
  crpx_tree_free (species);
  crpx_tree_free (gene);
}

static void
check_dates_and_shapes (void)
{
  crpx_tree *species, *gene;
  crpx_dtl_species *prepared = NULL;
  OK (crpx_tree_from_newick ("(A,B,C);", &species, NULL));
  assert (crpx_dtl_species_new (species, NULL, 0, &prepared) == CRPX_ERROR_UNSUPPORTED && !prepared);
  crpx_tree_free (species);
  OK (crpx_tree_from_newick ("((A,B),C);", &species, NULL));
  size_t nodes = crpx_tree_node_count (species);
  double ages[5];
  for (size_t node = 0; node < nodes; node++) {
    ages[node] = crpx_tree_first_child (species, (crpx_node_id) node) == CRPX_NODE_NONE ? 0 : 1;
    if (node == crpx_tree_root (species)) ages[node] = 2;
  }
  OK (crpx_dtl_species_new (species, ages, 3, &prepared));
  assert (crpx_dtl_species_states (prepared) == 6); /* C crosses the age-1 boundary. */
  crpx_dtl_workspace *workspace;
  OK (crpx_dtl_workspace_new (&workspace));
  crpx_dtl_costs costs = {10, 1, 1};
  OK (crpx_tree_from_newick ("((A|1,C|1),B|1);", &gene, NULL));
  crpx_node_id map[5];
  map_tips (gene, species, map);
  double value;
  OK (crpx_dtl_evaluate (gene, prepared, map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &value));
  crpx_dtl_mapping mapping[5];
  crpx_dtl_counts counts;
  OK (crpx_dtl_traceback (workspace, mapping, 5, &counts));
  assert (counts.transfers >= 1 && value == counts.transfers * costs.transfer + counts.losses * costs.loss +
          counts.duplications * costs.duplication);
  crpx_tree_free (gene);
  crpx_dtl_species_free (prepared);
  ages[crpx_tree_root (species)] = 1;
  assert (crpx_dtl_species_new (species, ages, 3, &prepared) == CRPX_ERROR_INVALID_ARGUMENT && !prepared);
  crpx_dtl_workspace_free (workspace);
  crpx_tree_free (species);
}

static void
check_temporal_restriction (void)
{
  crpx_tree *species, *gene;
  OK (crpx_tree_from_newick ("(((A,B),C),((D,E),F));", &species, NULL));
  OK (crpx_tree_from_newick ("(((B,A),F),(E,(D,C)));", &gene, NULL));
  /* AB and DE incoming branches occupy disjoint intervals. */
  double ages[11] = {0};
  crpx_node_id root = crpx_tree_root (species);
  crpx_node_id abc = crpx_tree_first_child (species, root);
  crpx_node_id def = crpx_tree_next_sibling (species, abc);
  ages[root] = 5;
  ages[abc] = 4;
  ages[crpx_tree_first_child (species, abc)] = 3;
  ages[def] = 2;
  ages[crpx_tree_first_child (species, def)] = 1;
  crpx_node_id tip_map[11];
  map_tips (gene, species, tip_map);
  crpx_dtl_species *dated, *undated;
  OK (crpx_dtl_species_new (species, ages, 6, &dated));
  OK (crpx_dtl_species_new (species, NULL, 0, &undated));
  crpx_dtl_workspace *workspace;
  OK (crpx_dtl_workspace_new (&workspace));
  crpx_dtl_costs costs = {10, 1, 1};
  double undated_cost, dated_cost;
  OK (crpx_dtl_evaluate (gene, undated, tip_map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &undated_cost));
  OK (crpx_dtl_evaluate (gene, dated, tip_map, &costs, CRPX_DTL_MINIMUM_COST, workspace, &dated_cost));
  assert (undated_cost == 3 && dated_cost == 4);
  crpx_dtl_mapping mapping[11];
  crpx_dtl_counts counts;
  OK (crpx_dtl_traceback (workspace, mapping, 11, &counts));
  assert (dated_cost == counts.duplications * costs.duplication + counts.transfers * costs.transfer +
          counts.losses * costs.loss);
  crpx_dtl_workspace_free (workspace);
  crpx_dtl_species_free (dated);
  crpx_dtl_species_free (undated);
  crpx_tree_free (species);
  crpx_tree_free (gene);
}

int
main (void)
{
  check_two_species (false);
  check_two_species (true);
  check_paralogs_and_losses ();
  check_dates_and_shapes ();
  check_temporal_restriction ();
  puts ("finite-history DTL: analytic partitions, event costs, paralogs, dates and workspace reuse passed");
  return 0;
}
