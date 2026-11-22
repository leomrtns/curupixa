/* SPDX-License-Identifier: MIT */
#include "inference.h"
#include "particle_weights.h"
#include "prob_distribution.h"
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
  crpx_tree_workspace *workspace;
  OK (crpx_tree_from_newick ("((A:1,B:1):2,C:3);", &tree, NULL));
  OK (crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace));
  double ages[5], rates[5], lengths[5], density;
  OK (crpx_tree_ages (tree, NULL, 1e-12, workspace, ages, 5));
  assert (ages[0] == 3 && ages[1] == 1 && ages[2] == 0);
  OK (crpx_clock_lengths (tree, ages, NULL, .5, lengths, 5));
  assert (lengths[1] == 1 && lengths[2] == .5 && lengths[4] == 1.5);
  double saved = lengths[1];
  ages[1] = 4;
  assert (crpx_clock_lengths (tree, ages, NULL, .5, lengths, 5) == CRPX_ERROR_INVALID_ARGUMENT && lengths[1] == saved);
  ages[1] = 1;
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 9876, 0);
  for (unsigned kind = 0; kind < 4; kind++) {
    crpx_clock_model model = {(crpx_clock_kind) kind, 2, .4};
    OK (crpx_clock_sample (&model, tree, ages, uniform, &stream, workspace, rates, 5));
    OK (crpx_clock_log_density (&model, tree, ages, rates, workspace, 5, &density));
    double expected = 0;
    for (size_t i = 1; i < 5; i++) {
      crpx_node_id parent = crpx_tree_parent (tree, (crpx_node_id) i);
      double mean = kind == 3 ? rates[parent] : 2;
      double v = .16;
      if (kind == 3) {
        crpx_node_id ancestor = crpx_tree_parent (tree, parent);
        v *= ((ages[parent] - ages[i]) + (ancestor == CRPX_NODE_NONE ? 0 : ages[ancestor] - ages[parent])) / 2;
      }
      if (kind == 1 || kind == 3) {
        double z = log (rates[i] / mean) + v / 2;
        expected += -log (rates[i]) - .5 * log (2 * acos (-1.0) * v) - z * z / (2 * v);
      } else if (kind == 2)
        expected += -log (2.0) - rates[i] / 2;
      else
        assert (rates[i] == 2);
    }
    assert (fabs (expected - density) < 1e-11);
  }
  /* Draw/scoring use the arithmetic mean, including autocorrelated expectations. */
  for (unsigned kind = 1; kind < 4; kind++) {
    crpx_clock_model model = {(crpx_clock_kind) kind, 2, .4};
    double sum = 0;
    for (size_t draw = 0; draw < 10000; draw++) {
      OK (crpx_clock_sample (&model, tree, ages, uniform, &stream, workspace, rates, 5));
      sum += rates[2];
    }
    assert (fabs (sum / 10000 - 2) < .07);
  }
  /* Heterochronous tip ages, zero-duration point masses and invalid inputs. */
  double tips[] = {0, 0, 0, 0, 2};
  assert (crpx_tree_ages (tree, tips, 1e-12, workspace, ages, 5) == CRPX_ERROR_INVALID_ARGUMENT);
  OK (crpx_tree_set_branch_length (tree, 4, 1));
  OK (crpx_tree_ages (tree, tips, 1e-12, workspace, ages, 5));
  assert (ages[0] == 3 && ages[4] == 2);
  OK (crpx_clock_lengths (tree, ages, NULL, .5, lengths, 5)); assert (lengths[4] == .5);
  OK (crpx_tree_set_branch_length (tree, 4, 3));
  OK (crpx_tree_ages (tree, NULL, 1e-12, workspace, ages, 5));
  crpx_clock_model strict = {CRPX_CLOCK_STRICT, 2, 0};
  OK (crpx_clock_sample (&strict, tree, ages, NULL, NULL, workspace, rates, 5));
  OK (crpx_clock_sample (&strict, tree, ages, uniform, &stream, workspace, rates, 5));
  rates[2] = 3;
  OK (crpx_clock_log_density (&strict, tree, ages, rates, workspace, 5, &density));
  assert (density == -INFINITY);
  crpx_tree_workspace_free (workspace);
  crpx_tree_free (tree);
  double logs[] = {-1000, -1000, -INFINITY}, weights[3], logsum, ess;
  OK (crpx_particle_normalise (logs, 3, weights, &logsum, &ess));
  assert (weights[0] == .5 && weights[2] == 0 && ess == 2 && fabs (logsum + 1000 - log (2.0)) < 1e-12);
  double scores[] = {-2, -1, 0}, equal[] = {1.0 / 3, 1.0 / 3, 1.0 / 3}, beta;
  OK (crpx_particle_temperature (equal, scores, 3, 0, .99, &beta));
  assert (beta > 0 && beta < 1);
  long double a = 0, b = 0;
  for (size_t i = 0; i < 3; i++) {
    double x = exp (beta * scores[i]);
    a += equal[i] * x;
    b += equal[i] * x * x;
  }
  assert (fabsl (a * a / b - .99) < 1e-12);
  size_t offspring[30];
  double sparse[] = {0, 1, 0};
  OK (crpx_particle_resample (sparse, 3, 30, false, uniform, &stream, offspring));
  for (size_t i = 0; i < 30; i++)
    assert (offspring[i] == 1);
  double zero_weights[] = {1, 0}, extreme[] = {-1e308, 1e308};
  OK (crpx_particle_temperature (zero_weights, extreme, 2, 0, .9, &beta));
  assert (beta == 1);
  return 0;
}
