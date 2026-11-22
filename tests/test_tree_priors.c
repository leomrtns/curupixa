/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "inference.h"
#include "global/random_number.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#define OK(call) assert ((call) == CRPX_OK)
static double
uniform (void *state)
{
  return crpx_rng_stream_double (state);
}
int
main (void)
{
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 42, 0);
  for (unsigned kind = 0; kind < 4; kind++)
    for (size_t leaves = 1; leaves < 16; leaves++) {
      crpx_tree_prior prior = {(crpx_tree_prior_kind) kind, 2, 1, 5, 100};
      crpx_tree *tree;
      OK (crpx_tree_prior_sample (&prior, leaves, NULL, uniform, &stream, &tree));
      assert (crpx_tree_leaf_count (tree) == leaves);
      OK (crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_BINARY));
      double value;
      OK (crpx_tree_prior_log_density (&prior, tree, &value));
      assert (isfinite (value));
      crpx_tree_free (tree);
    }
  crpx_tree *star, *resolved;
  OK (crpx_tree_from_newick ("(A:1,B:1,C:1,D:1,E:1);", &star, NULL));
  for (size_t i = 0; i < 100; i++) {
    OK (crpx_tree_resolve_random (star, uniform, &stream, &resolved));
    OK (crpx_tree_validate_policy (resolved, CRPX_TREE_REQUIRE_BINARY));
    assert (crpx_tree_node_count (resolved) == 9);
    crpx_tree_free (resolved);
  }
  crpx_tree_free (star);
  /* n=2: coalescent density is Exp(1/scale), including its sole labelled history. */
  crpx_tree *tree;
  OK (crpx_tree_from_newick ("(A:2,B:2);", &tree, NULL));
  crpx_tree_prior prior = {CRPX_TREE_PRIOR_COALESCENT, 1, 0, 5, 10};
  double value;
  OK (crpx_tree_prior_log_density (&prior, tree, &value));
  assert (fabs (value + log (10) + .2) < 1e-12);
  prior.kind = CRPX_TREE_PRIOR_UNIFORM;
  OK (crpx_tree_prior_log_density (&prior, tree, &value));
  assert (value == 0);
  crpx_tree_free (tree);
  /* Distribution checks distinguish PDA from the pair-coalescence topology law. */
  for (unsigned kind = 0; kind < 4; kind++) {
    crpx_tree_prior model = {(crpx_tree_prior_kind) kind, 2, kind == 2 ? 3 : 0, 2, 1};
    size_t balanced = 0;
    for (size_t sample = 0; sample < 4000; sample++) {
      OK (crpx_tree_prior_sample (&model, 4, NULL, uniform, &stream, &tree));
      crpx_node_id root = crpx_tree_root (tree), first = crpx_tree_first_child (tree, root);
      crpx_node_id second = crpx_tree_next_sibling (tree, first);
      balanced +=
        crpx_tree_first_child (tree, first) != CRPX_NODE_NONE && crpx_tree_first_child (tree, second) != CRPX_NODE_NONE;
      crpx_tree_free (tree);
    }
    double expected = kind == 0 ? .2 : 1.0 / 3;
    assert (fabs ((double) balanced / 4000 - expected) < .035);
  }
  /* Numerically integrate n=2 branching-time densities, including critical and subcritical BD. */
  for (unsigned model = 0; model < 4; model++) {
    crpx_tree_prior distribution = {CRPX_TREE_PRIOR_BIRTH_DEATH, 2, (double) model, 2, 1};
    double integral = 0;
    size_t steps = 2000;
    for (size_t i = 0; i < steps; i++) {
      double age = ((double) i + .5) * 2 / steps;
      OK (crpx_tree_new (&tree));
      crpx_node_id root, a, b;
      OK (crpx_tree_add_node (tree, CRPX_NODE_NONE, &root));
      OK (crpx_tree_add_node (tree, root, &a));
      OK (crpx_tree_add_node (tree, root, &b));
      OK (crpx_tree_set_branch_length (tree, a, age));
      OK (crpx_tree_set_branch_length (tree, b, age));
      OK (crpx_tree_prior_log_density (&distribution, tree, &value));
      integral += exp (value) * 2 / steps;
      crpx_tree_free (tree);
    }
    assert (fabs (integral - 1) < 1e-5);
  }
  OK (crpx_tree_from_newick ("(A,B,C,D);", &star, NULL));
  crpx_node_id constraint_tips[] = {1, 2, 3, 4};
  OK (crpx_tree_resolution_sample (star, uniform, &stream, &resolved));
  crpx_node_id candidate_tips[4];
  for (size_t node = 0; node < crpx_tree_node_count (resolved); node++) {
    const char *label = crpx_tree_label (resolved, (crpx_node_id) node);
    if (label && label[0] >= 'A' && label[0] <= 'D') candidate_tips[label[0] - 'A'] = (crpx_node_id) node;
  }
  OK (crpx_tree_resolution_log_density (star, constraint_tips, resolved, candidate_tips, 4, &value));
  assert (fabs (value + log (15)) < 1e-12);
  crpx_tree_free (star);
  crpx_tree_free (resolved);
  /* Parent import accepts arbitrary indices and rejects disconnected cycles. */
  crpx_node_id parents[] = {2, 2, CRPX_NODE_NONE, 0};
  const char *labels[] = {"internal", "B", "root", "A"};
  double lengths[] = {2, 3, .5, 1};
  OK (crpx_tree_from_parents (parents, 4, labels, lengths, &tree));
  OK (crpx_tree_validate (tree));
  assert (crpx_tree_leaf_count (tree) == 2);
  assert (crpx_tree_node_count (tree) == 4);
  double stem;
  OK (crpx_tree_branch_length (tree, crpx_tree_root (tree), &stem));
  assert (stem == .5);
  crpx_tree_free (tree);
  crpx_node_id cyclic[] = {CRPX_NODE_NONE, 2, 1};
  assert (crpx_tree_from_parents (cyclic, 3, NULL, NULL, &tree) == CRPX_ERROR_TOPOLOGY);
  assert (!tree);
  /* Wide constraints exercise linear construction and preserve their tip count. */
  size_t wide = 4001;
  crpx_node_id *wide_parents = calloc (wide, sizeof (*wide_parents));
  assert (wide_parents);
  wide_parents[0] = CRPX_NODE_NONE;
  OK (crpx_tree_from_parents (wide_parents, wide, NULL, NULL, &tree));
  OK (crpx_tree_resolution_sample (tree, uniform, &stream, &resolved));
  assert (crpx_tree_leaf_count (resolved) == wide - 1);
  OK (crpx_tree_validate_policy (resolved, CRPX_TREE_REQUIRE_BINARY));
  crpx_tree_free (resolved);
  crpx_tree_free (tree);
  free (wide_parents);
  return 0;
}
