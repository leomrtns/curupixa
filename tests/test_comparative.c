/* SPDX-License-Identifier: MIT */
#include "phylogenetics.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define OK(call) assert ((call) == CRPX_OK)
static crpx_node_id
leaf (const crpx_tree *tree, const char *label)
{
  for (size_t i = 0; i < crpx_tree_node_count (tree); i++) {
    const char *name = crpx_tree_label (tree, (crpx_node_id) i);
    if (name && !strcmp (name, label)) return (crpx_node_id) i;
  }
  assert (0);
  return CRPX_NODE_NONE;
}
static void
map (const crpx_tree *tree, size_t count, crpx_node_id *tips)
{
  for (size_t i = 0; i < count; i++) {
    char name[] = {(char) ('A' + i), 0};
    tips[i] = leaf (tree, name);
  }
}
static void
compare_distances (void)
{
  crpx_tree *first, *second, *star;
  OK (crpx_tree_from_newick ("((A:1,B:1):1,C:2);", &first, NULL));
  OK (crpx_tree_from_newick ("((A:2,B:1):1,C:2);", &second, NULL));
  OK (crpx_tree_from_newick ("(A:1,B:1,C:1);", &star, NULL));
  crpx_node_id a[3], b[3], c[3];
  map (first, 3, a);
  map (second, 3, b);
  map (star, 3, c);
  crpx_tree_comparison *x, *y;
  OK (crpx_tree_comparison_new (first, a, 3, true, true, &x));
  OK (crpx_tree_comparison_new (second, b, 3, true, true, &y));
  crpx_tree_distance_result result;
  OK (crpx_tree_comparison_distance (x, y, &result));
  assert (!result.rf && result.weighted_rf == 1 && result.branch_score == 1);
  double value;
  OK (crpx_tree_comparison_kendall_colijn (x, y, 0, &value));
  assert (value == 0);
  OK (crpx_tree_comparison_kendall_colijn (x, y, 1, &value));
  assert (value == 1);
  OK (crpx_tree_comparison_path_difference (x, y, true, &value));
  assert (fabs (value - sqrt (2)) < 1e-12);
  crpx_tree_comparison_free (x);
  crpx_tree_comparison_free (y);
  crpx_triplet_counts counts;
  OK (crpx_tree_triplet_distance (first, a, second, b, 3, &counts));
  assert (counts.resolved_agreement == 1);
  OK (crpx_tree_triplet_distance (first, a, star, c, 3, &counts));
  assert (counts.first_resolved_only == 1);
  crpx_tree_free (first);
  crpx_tree_free (second);
  crpx_tree_free (star);
  /* Equivalent unrooted storage roots combine their incident lengths. */
  OK (crpx_tree_from_newick ("((A:1,B:1):2,(C:1,D:1):3);", &first, NULL));
  OK (crpx_tree_from_newick ("(A:1,B:1,(C:1,D:1):5);", &second, NULL));
  OK (crpx_tree_weighted_distances (first, second, false, &result));
  assert (!result.rf && result.weighted_rf == 0 && result.branch_score == 0);
  OK (crpx_tree_path_difference (first, second, false, &value));
  assert (value == 0);
  crpx_tree_free (second);
  OK (crpx_tree_from_newick ("((A:1,C:1):2,(B:1,D:1):3);", &second, NULL));
  OK (crpx_tree_weighted_distances (first, second, false, &result));
  assert (result.rf == 2 && result.weighted_rf == 10);
  assert (fabs (result.branch_score - sqrt (50)) < 1e-12);
  size_t rf;
  OK (crpx_tree_rf_distance (first, second, false, &rf));
  assert (rf == result.rf);
  crpx_tree_free (second);
  OK (crpx_tree_from_newick ("(A:1,B:1,C:1,E:1);", &second, NULL));
  assert (crpx_tree_weighted_distances (first, second, false, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free (first);
  crpx_tree_free (second);
}
static void
regression (void)
{
  double identity[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double design[] = {1, 0, 1, 1, 1, 2, 1, 3}, response[] = {1, 2, 2, 4};
  crpx_covariance *covariance;
  crpx_regression_workspace *workspace;
  crpx_regression_result *result;
  OK (crpx_covariance_new (identity, 4, &covariance));
  OK (crpx_regression_workspace_new (&workspace));
  OK (crpx_gaussian_regression (covariance, design, 2, response, false, workspace, &result));
  assert (fabs (result->coefficients[0] - .9) < 1e-12 && fabs (result->coefficients[1] - .9) < 1e-12);
  assert (fabs (result->residual_sum_squares - .7) < 1e-12 && fabs (result->variance - .175) < 1e-12);
  assert (fabs (result->coefficient_covariance[0] - .1225) < 1e-12);
  crpx_regression_result_free (result);
  crpx_covariance *correlated;
  double matrix[] = {2, 1, 0, 1, 2, 1, 0, 1, 2}, intercept[] = {1, 1, 1}, observed[] = {1, 2, 4};
  OK (crpx_covariance_new (matrix, 3, &correlated));
  OK (crpx_gaussian_regression (correlated, intercept, 1, observed, false, workspace, &result));
  assert (fabs (result->coefficients[0] - 2.5) < 1e-12);
  assert (fabs (result->residual_sum_squares - 2.5) < 1e-12);
  assert (fabs (result->coefficient_covariance[0] - 5.0 / 6) < 1e-12);
  crpx_regression_result_free (result);
  crpx_covariance_free (correlated);
  crpx_regression_workspace_free (workspace);
  crpx_covariance_free (covariance);
}
static void
gaussian_predictions (void)
{
  crpx_tree *tree;
  crpx_node_id tips[2];
  OK (crpx_tree_from_newick ("(A:1,B:1);", &tree, NULL));
  map (tree, 2, tips);
  crpx_trait_covariance_options options;
  crpx_trait_covariance_options_init (&options);
  options.kind = CRPX_TRAIT_OU_STATIONARY;
  options.alpha = 1;
  crpx_gaussian_prediction *prediction;
  OK (crpx_phylogenetic_prediction (tree, tips, 2, &options, &prediction));
  double response[] = {1, 3}, observed_mean[] = {0, 0}, target_mean[] = {0, 0, 0}, means[3], variances[3];
  OK (crpx_gaussian_predict (prediction, response, observed_mean, target_mean, means, variances));
  double gain = exp (-1) / (1 + exp (-2));
  assert (fabs (means[0] - 4 * gain) < 1e-12 && fabs (variances[0] - (.5 - gain * exp (-1))) < 1e-12);
  assert (fabs (means[tips[0]] - 1) < 1e-12 && fabs (means[tips[1]] - 3) < 1e-12);
  crpx_gaussian_prediction_free (prediction);
  crpx_tree_free (tree);
  double identity[16] = {.1, 0, 0, 0, 0, .1, 0, 0, 0, 0, .1, 0, 0, 0, 0, .1};
  double design[] = {1, 1, 1, 1};
  uint8_t binary[] = {0, 0, 1, 1};
  crpx_covariance *covariance;
  crpx_regression_workspace *workspace;
  crpx_logistic_result *result;
  OK (crpx_covariance_new (identity, 4, &covariance));
  OK (crpx_regression_workspace_new (&workspace));
  OK (crpx_logistic_mixed_fit (covariance, design, 1, binary, 1e-10, 100, workspace, &result));
  assert (fabs (result->coefficients[0]) < 1e-10);
  for (size_t i = 0; i < 4; i++)
    assert (fabs (binary[i] - result->probabilities[i] - result->latent[i] / .1) < 1e-8);
  crpx_logistic_result_free (result);
  crpx_regression_workspace_free (workspace);
  crpx_covariance_free (covariance);
}
static void
quasi (void)
{
  crpx_quasi_sequence *sequence;
  double point[8];
  OK (crpx_quasi_sequence_new (CRPX_QUASI_SOBOL, 2, &sequence));
  const double expected[4][2] = {{.5, .5}, {.75, .25}, {.25, .75}, {.375, .375}};
  for (size_t i = 0; i < 4; i++) {
    OK (crpx_quasi_sequence_next (sequence, point, 8));
    assert (point[0] == expected[i][0] && point[1] == expected[i][1]);
  }
  OK (crpx_quasi_sequence_seek (sequence, 2));
  OK (crpx_quasi_sequence_next (sequence, point, 8));
  assert (point[0] == .75);
  crpx_quasi_sequence_free (sequence);
  OK (crpx_quasi_sequence_new (CRPX_QUASI_HALTON, 2, &sequence));
  OK (crpx_quasi_sequence_next (sequence, point, 8));
  assert (point[0] == .5 && fabs (point[1] - 1.0 / 3) < 1e-16);
  crpx_quasi_sequence_free (sequence);
  assert (crpx_quasi_sequence_new (CRPX_QUASI_SOBOL, 9, &sequence) == CRPX_ERROR_UNSUPPORTED && !sequence);
}
/* Brute-force all disjoint path pairings on four tips; no dynamic programming. */
static void
enumerate_pairs (const unsigned *states, uint64_t paths[4][4], unsigned remaining, uint64_t used, size_t pairs,
                 size_t support, crpx_pairing_result *best)
{
  if (!remaining) {
    if (pairs > best->pairs)
      *best = (crpx_pairing_result){pairs, support, support, 4};
    else if (pairs == best->pairs) {
      if (support < best->supporting_minimum) best->supporting_minimum = support;
      if (support > best->supporting_maximum) best->supporting_maximum = support;
    }
    return;
  }
  unsigned first = 0;
  while (!(remaining & (1u << first)))
    first++;
  remaining &= ~(1u << first);
  enumerate_pairs (states, paths, remaining, used, pairs, support, best);
  for (unsigned second = first + 1; second < 4; second++)
    if ((remaining & (1u << second)) && (states[first] ^ states[second]) == 3 && !(used & paths[first][second]))
      enumerate_pairs (states, paths, remaining & ~(1u << second), used | paths[first][second], pairs + 1,
                       support + (states[first] == 0 || states[first] == 3), best);
}
static void
association (void)
{
  crpx_tree *tree;
  crpx_node_id tips[4];
  crpx_association_workspace *workspace;
  OK (crpx_tree_from_newick ("((A:1,B:1):1,(C:1,D:1):1);", &tree, NULL));
  map (tree, 4, tips);
  uint64_t paths[4][4] = {{0}};
  for (size_t i = 0; i < 4; i++)
    for (size_t j = i + 1; j < 4; j++) {
      crpx_node_id ancestor;
      OK (crpx_tree_lca (tree, tips[i], tips[j], &ancestor));
      for (crpx_node_id x = tips[i]; x != ancestor; x = crpx_tree_parent (tree, x))
        paths[i][j] |= UINT64_C (1) << x;
      for (crpx_node_id x = tips[j]; x != ancestor; x = crpx_tree_parent (tree, x))
        paths[i][j] |= UINT64_C (1) << x;
    }
  OK (crpx_association_workspace_new (&workspace));
  for (unsigned code = 0; code < 256; code++) {
    uint8_t first[4], second[4];
    unsigned states[4];
    for (unsigned i = 0; i < 4; i++) {
      states[i] = (code >> (2 * i)) & 3;
      first[i] = (uint8_t) (states[i] >> 1);
      second[i] = (uint8_t) (states[i] & 1);
    }
    crpx_pairing_result result, oracle = {0, SIZE_MAX, 0, 4};
    enumerate_pairs (states, paths, 15, 0, 0, 0, &oracle);
    OK (crpx_trait_contrasting_pairs (tree, tips, 4, first, second, workspace, &result));
    assert (result.pairs == oracle.pairs && result.supporting_minimum == oracle.supporting_minimum &&
            result.supporting_maximum == oracle.supporting_maximum);
  }
  uint8_t clumped[] = {0, 0, 1, 1}, dispersed[] = {0, 1, 0, 1};
  double a, b, d;
  OK (crpx_trait_binary_dispersion (tree, tips, 4, clumped, workspace, &a));
  assert (a == 1);
  OK (crpx_trait_binary_dispersion (tree, tips, 4, dispersed, workspace, &b));
  assert (b == 2);
  OK (crpx_trait_d_statistic (a, a, b, &d));
  assert (d == 0);
  OK (crpx_trait_d_statistic (b, a, b, &d));
  assert (d == 1);
  crpx_association_workspace_free (workspace);
  crpx_tree_free (tree);
}
int
main (void)
{
  compare_distances ();
  regression ();
  gaussian_predictions ();
  quasi ();
  association ();
  puts ("comparative extensions passed");
  return 0;
}
