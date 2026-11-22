/* SPDX-License-Identifier: MIT */
#include "tanglegram.h"
#include <assert.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static crpx_node_id
tip (const crpx_tree *tree, const char *label)
{
  for (size_t node = 0; node < crpx_tree_node_count (tree); node++) {
    const char *name = crpx_tree_label (tree, (crpx_node_id) node);
    if (name && !strcmp (name, label)) return (crpx_node_id) node;
  }
  assert (0);
  return CRPX_NODE_NONE;
}

static uint64_t
count_crossings (const crpx_tanglegram_result *result, const crpx_tanglegram_link *links, size_t count)
{
  size_t first[64], second[64];
  uint64_t crossings = 0;
  for (size_t i = 0; i < result->first_leaves; i++) first[result->first_order[i]] = i;
  for (size_t i = 0; i < result->second_leaves; i++) second[result->second_order[i]] = i;
  for (size_t i = 0; i < count; i++) for (size_t j = 0; j < i; j++) {
    if (links[i].first == links[j].first || links[i].second == links[j].second) continue;
    crossings += (first[links[i].first] < first[links[j].first]) !=
                 (second[links[i].second] < second[links[j].second]);
  }
  return crossings;
}

int
main (void)
{
  crpx_tree *first, *second;
  crpx_tanglegram_workspace *workspace;
  crpx_tanglegram_result result;
  OK (crpx_tree_from_newick ("((A,B),(C,D));", &first, NULL));
  OK (crpx_tree_from_newick ("((A,C),(B,D));", &second, NULL));
  crpx_tanglegram_link links[5];
  const char *names[] = {"A", "B", "C", "D"};
  for (size_t i = 0; i < 4; i++) links[i] = (crpx_tanglegram_link) {tip (first, names[i]), tip (second, names[i])};
  uint64_t first_revision = crpx_tree_revision (first), second_revision = crpx_tree_revision (second);
  OK (crpx_tanglegram_workspace_new (first, second, links, 4, &workspace));
  OK (crpx_tanglegram_optimise (workspace, 1, &result));
  assert (result.evaluations == 1 && !result.exact);
  assert (count_crossings (&result, links, 4) == result.crossings);
  OK (crpx_tanglegram_optimise (workspace, 0, &result));
  assert (result.exact && result.crossings == 1);
  assert (count_crossings (&result, links, 4) == 1);
  assert (crpx_tree_revision (first) == first_revision && crpx_tree_revision (second) == second_revision);
  crpx_tanglegram_workspace_free (workspace);
  crpx_tree_free (second);

  /* Five tips versus three, including a polytomy and several copies per species. */
  crpx_tree_free (first);
  OK (crpx_tree_from_newick ("((A1,A2),(B,C,X));", &first, NULL));
  OK (crpx_tree_from_newick ("(C,B,A);", &second, NULL));
  links[0] = (crpx_tanglegram_link) {tip (first, "A1"), tip (second, "A")};
  links[1] = (crpx_tanglegram_link) {tip (first, "A2"), tip (second, "A")};
  links[2] = (crpx_tanglegram_link) {tip (first, "B"), tip (second, "B")};
  links[3] = (crpx_tanglegram_link) {tip (first, "C"), tip (second, "C")};
  links[4] = links[0];
  OK (crpx_tanglegram_workspace_new (first, second, links, 5, &workspace));
  crpx_tree_free (first);
  crpx_tree_free (second);
  OK (crpx_tanglegram_optimise (workspace, 1000, &result));
  assert (result.exact && result.crossings == 0 && result.first_leaves == 5 && result.second_leaves == 3);
  assert (count_crossings (&result, links, 5) == 0);
  crpx_tanglegram_workspace_free (workspace);

  /* Complete bipartite links force one crossing per pair of tips on each side. */
  OK (crpx_tree_from_newick ("(A,B);", &first, NULL));
  OK (crpx_tree_from_newick ("(X,Y,Z);", &second, NULL));
  crpx_tanglegram_link complete[6];
  for (size_t i = 0; i < 2; i++) for (size_t j = 0; j < 3; j++)
    complete[i * 3 + j] = (crpx_tanglegram_link) {(crpx_node_id) (i + 1), (crpx_node_id) (j + 1)};
  OK (crpx_tanglegram_workspace_new (first, second, complete, 6, &workspace));
  OK (crpx_tanglegram_optimise (workspace, 0, &result));
  assert (result.exact && result.crossings == 3);
  assert (count_crossings (&result, complete, 6) == 3);
  crpx_tanglegram_workspace_free (workspace);
  /* Independent enumeration of all 2! * 3! embeddings for every simple
   * bipartite link subset; compare both the initial counter and optimum. */
  const crpx_node_id right_orders[6][3] = {{1,2,3}, {1,3,2}, {2,1,3}, {2,3,1}, {3,1,2}, {3,2,1}};
  const crpx_node_id left_orders[2][2] = {{1,2}, {2,1}};
  for (unsigned mask = 0; mask < 64; mask++) {
    crpx_tanglegram_link subset[6];
    size_t count = 0;
    for (size_t i = 0; i < 6; i++) if (mask & (1u << i)) subset[count++] = complete[i];
    uint64_t minimum = UINT64_MAX;
    for (size_t i = 0; i < 2; i++) for (size_t j = 0; j < 6; j++) {
      crpx_tanglegram_result candidate = {0};
      candidate.first_order = left_orders[i]; candidate.first_leaves = 2;
      candidate.second_order = right_orders[j]; candidate.second_leaves = 3;
      uint64_t value = count_crossings (&candidate, subset, count);
      if (value < minimum) minimum = value;
    }
    OK (crpx_tanglegram_workspace_new (first, second, subset, count, &workspace));
    OK (crpx_tanglegram_optimise (workspace, 1, &result));
    assert (result.crossings == count_crossings (&result, subset, count));
    OK (crpx_tanglegram_optimise (workspace, 0, &result));
    assert (result.exact && result.crossings == minimum);
    crpx_tanglegram_workspace_free (workspace);
  }
  OK (crpx_tanglegram_workspace_new (first, second, NULL, 0, &workspace));
  OK (crpx_tanglegram_optimise (workspace, 1, &result));
  assert (result.exact && result.crossings == 0);
  crpx_tanglegram_workspace_free (workspace);
  complete[0].first = crpx_tree_root (first);
  assert (crpx_tanglegram_workspace_new (first, second, complete, 6, &workspace) == CRPX_ERROR_INVALID_ARGUMENT);
  assert (!workspace);
  crpx_tree_free (first);
  crpx_tree_free (second);
  return 0;
}
