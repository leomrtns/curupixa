/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "tree_search.h"
#include "topology_distance.h"
#include <assert.h>
#define OK(call) assert ((call) == CRPX_OK)
static crpx_status
rf_score (const crpx_tree *tree, void *target, double *score)
{
  size_t distance;
  crpx_status status = crpx_tree_rf_distance (tree, target, true, &distance);
  if (status == CRPX_OK) *score = (double) distance;
  return status;
}
static crpx_status
unrooted_rf_score (const crpx_tree *tree, void *target, double *score)
{
  size_t distance;
  crpx_status status = crpx_tree_rf_distance (tree, target, false, &distance);
  if (status == CRPX_OK) *score = (double) distance;
  return status;
}
int
main (void)
{
  crpx_tree *start, *target;
  crpx_tree_search_result result;
  OK (crpx_tree_from_newick ("((A,B),(C,D));", &start, NULL));
  OK (crpx_tree_from_newick ("(A,(B,(C,D)));", &target, NULL));
  uint64_t revision = crpx_tree_revision (start);
  OK (crpx_tree_search_nni (start, rf_score, target, 0, &result));
  assert (result.local_optimum && result.score == 0 && result.moves > 0);
  assert (crpx_tree_revision (start) == revision);
  crpx_tree_free (result.tree);
  OK (crpx_tree_search_nni (target, rf_score, target, 0, &result));
  assert (result.local_optimum && result.score == 0 && result.moves == 0);
  crpx_tree_free (result.tree);
  crpx_tree_free (target);
  OK (crpx_tree_from_newick ("((A,C),(B,D));", &target, NULL));
  OK (crpx_tree_search_unrooted_nni (start, unrooted_rf_score, target, 0, &result));
  assert (result.local_optimum && result.score == 0 && result.moves == 1);
  assert (crpx_tree_revision (start) == revision);
  crpx_tree_free (result.tree);
  crpx_tree_free (start);
  crpx_tree_free (target);
  return 0;
}
