/* SPDX-License-Identifier: MIT */
#include "reconstruction.h"
#include "stochastic_mapping.h"
#include "global/random_number.h"
#include <assert.h>
#include <math.h>
#define OK(call) assert ((call) == CRPX_OK)
static double
uniform (void *state)
{
  return crpx_rng_stream_double (state);
}
int
main (void)
{
  crpx_tree *tree;
  OK (crpx_tree_from_newick ("(A:.4,B:.7);", &tree, NULL));
  const char *names[] = {"A", "B"};
  double observations[] = {1, .2, .3, 0, .1, .4, 0, 1};
  crpx_character_data *data;
  OK (crpx_character_data_from_tip_likelihoods (names, 2, 1, 4, observations, NULL, false, &data));
  const double rates[] = {1, 2, 3, 1, .5, 1.5, 2, 1}, root[] = {.1, .2, .3, .4};
  crpx_substitution_model *model;
  OK (crpx_substitution_model_binary_pair (rates, root, &model));
  double left[16], right[16];
  OK (crpx_substitution_model_transition (model, .4, CRPX_BRANCH_FIXED, left, 16));
  OK (crpx_substitution_model_transition (model, .7, CRPX_BRANCH_FIXED, right, 16));
  double sum = 0, maximum = 0, marginal[3][4] = {{0}};
  for (size_t r = 0; r < 4; r++)
    for (size_t a = 0; a < 4; a++)
      for (size_t b = 0; b < 4; b++) {
        double value = root[r] * left[4 * r + a] * right[4 * r + b] * observations[a] * observations[4 + b];
        sum += value;
        if (value > maximum) maximum = value;
        marginal[0][r] += value;
        marginal[1][a] += value;
        marginal[2][b] += value;
      }
  crpx_reconstruction_workspace *workspace;
  OK (crpx_reconstruction_workspace_new (&workspace));
  crpx_likelihood_reconstruction *reconstruction;
  OK (crpx_likelihood_reconstruct (tree, data, model, NULL, workspace, &reconstruction));
  assert (fabs (crpx_likelihood_reconstruction_score (reconstruction) - log (sum)) < 1e-11);
  for (size_t node = 0; node < 3; node++)
    for (size_t state = 0; state < 4; state++)
      assert (fabs (crpx_likelihood_reconstruction_probabilities (reconstruction, (crpx_node_id) node, 0)[state] -
                    marginal[node][state] / sum) < 1e-11);
  crpx_likelihood_reconstruction_free (reconstruction);
  crpx_joint_reconstruction *joint;
  OK (crpx_likelihood_joint_states (tree, data, model, NULL, false, NULL, NULL, workspace, &joint));
  assert (fabs (joint->pattern_log_joint[0] - log (maximum)) < 1e-11);
  crpx_joint_reconstruction_free (joint);
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 123, 0);
  size_t counts[4] = {0};
  for (size_t i = 0; i < 2000; i++) {
    OK (crpx_likelihood_joint_states (tree, data, model, NULL, true, uniform, &stream, workspace, &joint));
    counts[joint->assignment[0]]++;
    crpx_joint_reconstruction_free (joint);
  }
  for (size_t state = 0; state < 4; state++)
    assert (fabs ((double) counts[state] / 2000 - marginal[0][state] / sum) < .05);
  crpx_reconstruction_workspace_free (workspace);
  crpx_character_data_free (data);
  crpx_tree_free (tree);
  crpx_mapping_workspace *mapping;
  OK (crpx_mapping_workspace_new (&mapping));
  crpx_branch_history *history;
  OK (crpx_substitution_bridge_sample (model, .5, 0, 3, 1e-10, 200, uniform, &stream, mapping, &history));
  assert (history->count >= 2 && history->states[0] == 0 && history->states[history->count] == 3);
  for (size_t i = 0; i < history->count; i++) {
    size_t change = history->states[i] ^ history->states[i + 1];
    assert (change == 1 || change == 2);
    assert (history->times[i] > 0 && history->times[i] < .5);
    if (i) assert (history->times[i] >= history->times[i - 1]);
  }
  crpx_branch_history_free (history);
  crpx_substitution_model_free (model);
  OK (crpx_substitution_model_poisson (2, &model));
  size_t total = 0;
  for (size_t i = 0; i < 2000; i++) {
    OK (crpx_substitution_bridge_sample (model, 1, 0, 1, 1e-10, 100, uniform, &stream, mapping, &history));
    assert (history->count % 2 == 1);
    total += history->count;
    crpx_branch_history_free (history);
  }
  assert (fabs ((double) total / 2000 - 1 / tanh (1)) < .08);
  assert (crpx_substitution_bridge_sample (model, 100, 0, 1, 1e-10, 1, uniform, &stream, mapping, &history) ==
          CRPX_ERROR_CONVERGENCE);
  crpx_mapping_workspace_free (mapping);
  crpx_substitution_model_free (model);
  const double independent_rates[] = {1, 1, 1, 1, 1, 1, 1, 1}, uniform_root[] = {1, 1, 1, 1};
  OK (crpx_substitution_model_binary_pair (independent_rates, uniform_root, &model));
  double matrix[16];
  OK (crpx_substitution_model_transition (model, .7, CRPX_BRANCH_FIXED, matrix, 16));
  double flip = .5 * (1 - exp (-1.4));
  assert (fabs (matrix[0] - (1 - flip) * (1 - flip)) < 1e-13 && fabs (matrix[3] - flip * flip) < 1e-13);
  OK (crpx_substitution_model_transition (model, .7, CRPX_BRANCH_SITEWISE_EXPONENTIAL, matrix, 16));
  assert (fabs (matrix[3] - .25 * (1 - 2 / 2.4 + 1 / 3.8)) < 1e-13);
  OK (crpx_substitution_model_transition (model, 1e200, CRPX_BRANCH_FIXED, matrix, 16));
  for (size_t i = 0; i < 16; i++)
    assert (fabs (matrix[i] - .25) < 1e-12);
  crpx_substitution_model_free (model);
  return 0;
}
