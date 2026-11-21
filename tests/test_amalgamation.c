/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "amalgamation.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)

static void
map_tips (const crpx_tree *tree, const char *const *names, size_t taxa, crpx_node_id *tips)
{
  for (size_t i = 0; i < taxa; i++) {
    tips[i] = CRPX_NODE_NONE;
    for (size_t j = 0; j < crpx_tree_node_count (tree); j++) {
      const char *label = crpx_tree_label (tree, (crpx_node_id) j);
      if (label && !strcmp (label, names[i])) tips[i] = (crpx_node_id) j;
    }
    assert (tips[i] != CRPX_NODE_NONE);
  }
}

static double
uniform_draw (void *state)
{
  uint64_t *seed = state;
  *seed = *seed * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
  return (double) (*seed >> 11) * 0x1p-53;
}

static double bad_draw (void *state) { (void) state; return 1; }

static void
check_distribution (void)
{
  const char *names[] = {"A", "B", "C", "D", "E", "F"};
  const char *newick[] = {
    "(((A,B),C),((D,E),F));", "(((A,C),B),((D,F),E));",
    "(((A,B),C),((D,F),E));", "(((A,C),B),((D,E),F));",
    "(((B,C),A),((D,E),F));"
  };
  crpx_tree *trees[5];
  crpx_node_id tips[5][6];
  for (size_t i = 0; i < 5; i++) {
    OK (crpx_tree_from_newick (newick[i], trees + i, NULL));
    map_tips (trees[i], names, 6, tips[i]);
  }
  const crpx_tree *sample[] = {trees[0], trees[1]};
  const crpx_node_id *sample_tips[] = {tips[0], tips[1]};
  crpx_amalgamation *model;
  crpx_amalgamation_workspace *workspace;
  OK (crpx_amalgamation_new (sample, sample_tips, 2, 6, names, NULL, &model));
  OK (crpx_amalgamation_workspace_new (&workspace));
  assert (crpx_amalgamation_taxa (model) == 6 && crpx_amalgamation_samples (model) == 2);
  double sum = 0;
  for (size_t i = 0; i < 4; i++) {
    double probability;
    OK (crpx_amalgamation_log_probability (model, trees[i], tips[i], workspace, &probability));
    assert (fabs (probability - log (0.25)) < 1e-12);
    sum += exp (probability);
  }
  assert (fabs (sum - 1) < 1e-12);
  double probability;
  OK (crpx_amalgamation_log_probability (model, trees[4], tips[4], workspace, &probability));
  assert (probability == -INFINITY);
  uint64_t seed = 1729;
  size_t unseen_samples = 0;
  for (size_t i = 0; i < 100; i++) {
    crpx_tree *generated;
    crpx_node_id generated_tips[6];
    double expected;
    OK (crpx_amalgamation_sample (model, uniform_draw, &seed, &generated, &probability));
    map_tips (generated, names, 6, generated_tips);
    OK (crpx_amalgamation_log_probability (model, generated, generated_tips, workspace, &expected));
    assert (fabs (expected - probability) < 1e-12 && fabs (expected - log (0.25)) < 1e-12);
    /* C and B are root children of the left clade for different input trees;
     * checking the two cherry pairs detects the new combinations. */
    bool ab = crpx_tree_parent (generated, generated_tips[0]) == crpx_tree_parent (generated, generated_tips[1]);
    bool de = crpx_tree_parent (generated, generated_tips[3]) == crpx_tree_parent (generated, generated_tips[4]);
    unseen_samples += ab != de;
    crpx_tree_free (generated);
  }
  assert (unseen_samples > 20);
  crpx_tree *generated = (void *) 1;
  assert (crpx_amalgamation_sample (model, bad_draw, NULL, &generated, NULL) == CRPX_ERROR_INVALID_ARGUMENT && !generated);
  crpx_amalgamation_free (model);
  uint64_t weights[] = {3, 1};
  OK (crpx_amalgamation_new (sample, sample_tips, 2, 6, names, weights, &model));
  OK (crpx_amalgamation_log_probability (model, trees[0], tips[0], workspace, &probability));
  assert (fabs (probability - log (0.75 * 0.75)) < 1e-12);
  crpx_amalgamation_free (model);
  crpx_tree *star;
  OK (crpx_tree_from_newick ("(A,B,C,D,E,F);", &star, NULL));
  crpx_node_id star_tips[6]; map_tips (star, names, 6, star_tips);
  const crpx_tree *bad_sample[] = {star}; const crpx_node_id *bad_tips[] = {star_tips};
  assert (crpx_amalgamation_new (bad_sample, bad_tips, 1, 6, names, NULL, &model) == CRPX_ERROR_UNSUPPORTED && !model);
  crpx_tree_free (star);
  crpx_amalgamation_workspace_free (workspace);
  for (size_t i = 0; i < 5; i++) crpx_tree_free (trees[i]);
}

static void
check_large_clade (void)
{
  enum {taxa = 65};
  crpx_tree *tree, *generated;
  crpx_node_id tips[taxa], parent, next;
  char storage[taxa][12]; const char *names[taxa];
  OK (crpx_tree_new (&tree));
  OK (crpx_tree_add_node (tree, CRPX_NODE_NONE, &parent));
  for (size_t i = 0; i < taxa; i++) {
    snprintf (storage[i], sizeof (storage[i]), "copy%zu", i); names[i] = storage[i];
    OK (crpx_tree_add_node (tree, parent, tips + i));
    OK (crpx_tree_set_label (tree, tips[i], names[i]));
    if (i + 2 < taxa) { OK (crpx_tree_add_node (tree, parent, &next)); parent = next; }
  }
  const crpx_tree *sample[] = {tree}; const crpx_node_id *maps[] = {tips};
  crpx_amalgamation *model;
  OK (crpx_amalgamation_new (sample, maps, 1, taxa, names, NULL, &model));
  uint64_t seed = 3; double probability;
  OK (crpx_amalgamation_sample (model, uniform_draw, &seed, &generated, &probability));
  assert (probability == 0 && crpx_tree_leaf_count (generated) == taxa);
  crpx_tree_free (generated); crpx_tree_free (tree); crpx_amalgamation_free (model);
}

int main (void)
{
  check_distribution ();
  check_large_clade ();
  puts ("amalgamation tests passed");
  return 0;
}
