/* Full-data messages inspired by branchtraits/src/discrete.c, generalised to
 * curupixa model/data contracts. Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Reconstruct marginal ancestral probabilities or optimal parsimony states from discrete observations. */
#include "reconstruction.h"
#include "parsimony_internal.h"
#include "likelihood_kernel_internal.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_reconstruction_workspace {
  crpx_likelihood_workspace *likelihood;
  crpx_tree_workspace *traversal;
  double *buffers;
  size_t *zeros;
  crpx_node_id *tips;
  size_t buffer_capacity, zero_capacity, tip_capacity;
};
struct crpx_likelihood_reconstruction {
  size_t nodes, patterns, states;
  double score, *probabilities;
};
struct crpx_parsimony_reconstruction {
  size_t nodes, patterns, states;
  double score, *pattern_scores;
  size_t *assignment;
  unsigned char *optimal;
};
crpx_status
crpx_reconstruction_workspace_new (crpx_reconstruction_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void
crpx_reconstruction_workspace_clear (crpx_reconstruction_workspace *w)
{
  if (!w) return;
  crpx_likelihood_workspace_free (w->likelihood);
  crpx_tree_workspace_free (w->traversal);
  free (w->buffers);
  free (w->zeros);
  free (w->tips);
  memset (w, 0, sizeof (*w));
}
void
crpx_reconstruction_workspace_free (crpx_reconstruction_workspace *w)
{
  if (w) {
    crpx_reconstruction_workspace_clear (w);
    free (w);
  }
}
static void *
grow (void *array, size_t *capacity, size_t count, size_t unit, crpx_status *status)
{
  if (*status != CRPX_OK || count <= *capacity) return array;
  if (count > SIZE_MAX / unit) {
    *status = CRPX_ERROR_OVERFLOW;
    return array;
  }
  void *p = realloc (array, count * unit);
  if (!p) {
    *status = CRPX_ERROR_ALLOCATION;
    return array;
  }
  *capacity = count;
  return p;
}
static crpx_status
prepare (const crpx_tree *tree, const crpx_character_data *data, crpx_reconstruction_workspace *w, bool kernels)
{
  size_t n = crpx_tree_node_count (tree), k = crpx_character_data_states (data), t = crpx_character_data_taxa (data);
  if (!n || !k || !t) return CRPX_ERROR_INVALID_ARGUMENT;
  if (k > SIZE_MAX / 5 || n > SIZE_MAX / k || n * k > (SIZE_MAX - 5 * k) / 3) return CRPX_ERROR_OVERFLOW;
  crpx_status status = CRPX_OK;
  w->buffers = grow (w->buffers, &w->buffer_capacity, 3 * n * k + 5 * k, sizeof (double), &status);
  w->zeros = grow (w->zeros, &w->zero_capacity, k, sizeof (size_t), &status);
  if (!kernels) w->tips = grow (w->tips, &w->tip_capacity, t, sizeof (crpx_node_id), &status);
  if (status == CRPX_OK && !w->traversal)
    status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &w->traversal);
  /* Likelihood evaluation already validates/maps observations. Parsimony owns that setup here; avoid a redundant map in
   the likelihood reconstruction. */
  if (status == CRPX_OK && !kernels) status = crpx_character_data_map_tips (data, tree, w->tips, t);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (w->traversal, tree);
  return status;
}
static void
normalise (double *v, size_t k)
{
  double maximum = -INFINITY;
  for (size_t s = 0; s < k; ++s)
    if (v[s] > maximum) maximum = v[s];
  if (isfinite (maximum))
    for (size_t s = 0; s < k; ++s)
      v[s] -= maximum;
}
static double
log_sum (const double *v, size_t k)
{
  double value = -INFINITY;
  for (size_t s = 0; s < k; ++s)
    value = crpx_private_likelihood_log_add (value, v[s]);
  return value;
}
void
crpx_likelihood_reconstruction_free (crpx_likelihood_reconstruction *r)
{
  if (r) {
    free (r->probabilities);
    free (r);
  }
}
size_t
crpx_likelihood_reconstruction_nodes (const crpx_likelihood_reconstruction *r)
{
  return r ? r->nodes : 0;
}
size_t
crpx_likelihood_reconstruction_patterns (const crpx_likelihood_reconstruction *r)
{
  return r ? r->patterns : 0;
}
size_t
crpx_likelihood_reconstruction_states (const crpx_likelihood_reconstruction *r)
{
  return r ? r->states : 0;
}
double
crpx_likelihood_reconstruction_score (const crpx_likelihood_reconstruction *r)
{
  return r ? r->score : NAN;
}
const double *
crpx_likelihood_reconstruction_probabilities (const crpx_likelihood_reconstruction *r, crpx_node_id node,
                                              size_t pattern)
{
  return r && node < r->nodes && pattern < r->patterns ? r->probabilities + (pattern * r->nodes + node) * r->states
                                                       : NULL;
}
crpx_status
crpx_likelihood_reconstruct (const crpx_tree *tree, const crpx_character_data *data, const crpx_substitution_model *m,
                             const crpx_likelihood_options *options, crpx_reconstruction_workspace *w,
                             crpx_likelihood_reconstruction **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !data || !m || !w) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = crpx_tree_node_count (tree), k = crpx_substitution_model_states (m),
         p = crpx_character_data_patterns (data);
  if (!n || !k || k != crpx_character_data_states (data)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n > SIZE_MAX / k || (p && n * k > SIZE_MAX / p / sizeof (double))) return CRPX_ERROR_OVERFLOW;
  crpx_status status = prepare (tree, data, w, true);
  crpx_likelihood_reconstruction *r = NULL;
  if (status != CRPX_OK) goto done;
  if (!w->likelihood) status = crpx_likelihood_workspace_new (&w->likelihood);
  if (status != CRPX_OK) goto done;
  crpx_likelihood_options o;
  crpx_likelihood_options_init (&o);
  if (options) o = *options;
  o.retain_partials = true;
  double score;
  status = crpx_likelihood_evaluate (tree, data, m, &o, w->likelihood, &score);
  if (status != CRPX_OK) goto done;
  r = calloc (1, sizeof (*r));
  if (!r) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  r->nodes = n;
  r->states = k;
  r->patterns = p;
  r->score = score;
  r->probabilities = p ? malloc (p * n * k * sizeof (double)) : NULL;
  if (p && !r->probabilities) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  double *down = w->buffers, *up = down + n * k, *message = up + n * k, *scratch = message + n * k;
  double *total = scratch + k, *exclude = total + k, *temp = exclude + k, *log_pi = temp + k;
  const double *pi = crpx_substitution_model_frequencies (m);
  for (size_t s = 0; s < k; ++s)
    log_pi[s] = log (pi[s]);
  crpx_node_id root = crpx_tree_root (tree);
  const crpx_log_kernel *kernels = crpx_private_likelihood_kernels (w->likelihood);
  const crpx_node_id *pre = crpx_tree_workspace_preorder (w->traversal);
  const double *pattern_logs = crpx_likelihood_workspace_pattern_log_likelihoods (w->likelihood, NULL);
  for (size_t pattern = 0; pattern < p; ++pattern) {
    if (!isfinite (pattern_logs[pattern])) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    for (size_t node = 0; node < n; ++node) {
      const double *partial;
      double scale;
      status = crpx_likelihood_workspace_log_partial (w->likelihood, (crpx_node_id) node, pattern, &partial, &scale);
      if (status != CRPX_OK) goto done;
      memcpy (down + node * k, partial, k * sizeof (double));
      if (node != root) {
        crpx_private_likelihood_project (m, kernels + node, partial, log_pi, message + node * k);
        normalise (message + node * k, k);
      }
    }
    memcpy (up + (size_t) root * k, log_pi, k * sizeof (double));
    for (size_t pos = 0; pos < n; ++pos) {
      crpx_node_id node = pre[pos];
      for (size_t s = 0; s < k; ++s) {
        total[s] = 0;
        w->zeros[s] = 0;
      }
      for (crpx_node_id child = crpx_tree_first_child (tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling (tree, child))
        for (size_t s = 0; s < k; ++s) {
          double value = message[(size_t) child * k + s];
          if (value == -INFINITY)
            ++w->zeros[s];
          else
            total[s] += value;
        }
      for (crpx_node_id child = crpx_tree_first_child (tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling (tree, child)) {
        for (size_t s = 0; s < k; ++s) {
          double own = message[(size_t) child * k + s];
          size_t zeros = w->zeros[s] - (own == -INFINITY);
          exclude[s] = zeros ? -INFINITY : up[(size_t) node * k + s] + total[s] - (isfinite (own) ? own : 0);
        }
        crpx_private_likelihood_project_transpose (m, kernels + child, exclude, log_pi, temp, up + (size_t) child * k);
        normalise (up + (size_t) child * k, k);
      }
      for (size_t s = 0; s < k; ++s)
        scratch[s] = down[(size_t) node * k + s] + up[(size_t) node * k + s];
      double norm = log_sum (scratch, k);
      if (!isfinite (norm)) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      double *prob = r->probabilities + (pattern * n + node) * k;
      for (size_t s = 0; s < k; ++s)
        prob[s] = exp (scratch[s] - norm);
    }
  }
  *out = r;
  r = NULL;
done:
  crpx_likelihood_reconstruction_free (r);
  crpx_tree_workspace_reset (w->traversal);
  return status;
}
void
crpx_parsimony_reconstruction_free (crpx_parsimony_reconstruction *r)
{
  if (r) {
    free (r->assignment);
    free (r->optimal);
    free (r->pattern_scores);
    free (r);
  }
}
size_t
crpx_parsimony_reconstruction_nodes (const crpx_parsimony_reconstruction *r)
{
  return r ? r->nodes : 0;
}
size_t
crpx_parsimony_reconstruction_patterns (const crpx_parsimony_reconstruction *r)
{
  return r ? r->patterns : 0;
}
size_t
crpx_parsimony_reconstruction_states (const crpx_parsimony_reconstruction *r)
{
  return r ? r->states : 0;
}
double
crpx_parsimony_reconstruction_score (const crpx_parsimony_reconstruction *r)
{
  return r ? r->score : NAN;
}
const double *
crpx_parsimony_reconstruction_pattern_scores (const crpx_parsimony_reconstruction *r)
{
  return r ? r->pattern_scores : NULL;
}
size_t
crpx_parsimony_reconstruction_assignment (const crpx_parsimony_reconstruction *r, crpx_node_id node, size_t p)
{
  return r && node < r->nodes && p < r->patterns ? r->assignment[p * r->nodes + node] : SIZE_MAX;
}
const unsigned char *
crpx_parsimony_reconstruction_optimal_states (const crpx_parsimony_reconstruction *r, crpx_node_id node, size_t pattern)
{
  return r && node < r->nodes && pattern < r->patterns ? r->optimal + (pattern * r->nodes + node) * r->states : NULL;
}
static double
change (size_t a, size_t b, crpx_parsimony_cost_mode mode, const double *costs, size_t k)
{
  if (mode == CRPX_PARSIMONY_CUSTOM) return costs[a * k + b];
  if (mode == CRPX_PARSIMONY_ORDERED) return (double) (a > b ? a - b : b - a);
  return a != b;
}
crpx_status
crpx_parsimony_reconstruct (const crpx_tree *tree, const crpx_character_data *data, crpx_parsimony_cost_mode mode,
                            const double *costs, size_t cost_count, double tolerance, crpx_reconstruction_workspace *w,
                            crpx_parsimony_reconstruction **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !data || !w || !isfinite (tolerance) || tolerance < 0 ||
      (mode != CRPX_PARSIMONY_UNORDERED && mode != CRPX_PARSIMONY_ORDERED && mode != CRPX_PARSIMONY_CUSTOM))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n = crpx_tree_node_count (tree), k = crpx_character_data_states (data),
         p = crpx_character_data_patterns (data);
  size_t taxa = crpx_character_data_taxa (data);
  if (!n || !k) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n > SIZE_MAX / k || (p && (n > SIZE_MAX / p / sizeof (size_t) || n * k > SIZE_MAX / p)) ||
      p > SIZE_MAX / sizeof (double))
    return CRPX_ERROR_OVERFLOW;
  if (mode == CRPX_PARSIMONY_CUSTOM) {
    if (k > SIZE_MAX / k) return CRPX_ERROR_OVERFLOW;
    if (!costs || cost_count != k * k) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t i = 0; i < cost_count; ++i)
      if (!isfinite (costs[i]) || costs[i] < 0 || (i / k == i % k && costs[i] != 0)) return CRPX_ERROR_INVALID_ARGUMENT;
  } else if (costs || cost_count)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare (tree, data, w, false);
  crpx_parsimony_reconstruction *r = NULL;
  if (status != CRPX_OK) goto done;
  r = calloc (1, sizeof (*r));
  if (!r) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  r->nodes = n;
  r->states = k;
  r->patterns = p;
  r->assignment = p ? malloc (p * n * sizeof (size_t)) : NULL;
  r->optimal = p ? malloc (p * n * k) : NULL;
  r->pattern_scores = p ? malloc (p * sizeof (double)) : NULL;
  if (p && (!r->assignment || !r->optimal || !r->pattern_scores)) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  double *down = w->buffers, *up = down + n * k, *message = up + n * k, *scratch = message + n * k,
         *total = scratch + k;
  const double *values = crpx_character_data_pattern_data (data);
  const size_t *weights = crpx_character_data_pattern_weights (data);
  const crpx_node_id *post = crpx_tree_workspace_postorder (w->traversal),
                     *pre = crpx_tree_workspace_preorder (w->traversal);
  crpx_node_id root = crpx_tree_root (tree);
  long double weighted = 0;
  for (size_t pattern = 0; pattern < p; ++pattern) {
    memset (down, 0, n * k * sizeof (double));
    memset (up, 0, n * k * sizeof (double));
    for (size_t t = 0; t < taxa; ++t)
      for (size_t s = 0; s < k; ++s)
        down[(size_t) w->tips[t] * k + s] = values[(pattern * taxa + t) * k + s] > 0 ? 0 : INFINITY;
    for (size_t pos = 0; pos < n; ++pos) {
      crpx_node_id node = post[pos];
      for (crpx_node_id child = crpx_tree_first_child (tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling (tree, child))
        for (size_t s = 0; s < k; ++s)
          down[(size_t) node * k + s] =
            crpx_private_parsimony_add (down[(size_t) node * k + s], message[(size_t) child * k + s]);
      crpx_private_parsimony_transform (down + (size_t) node * k, message + (size_t) node * k, k, mode, costs, false);
    }
    double optimum = INFINITY;
    for (size_t s = 0; s < k; ++s)
      if (down[(size_t) root * k + s] < optimum) {
        optimum = down[(size_t) root * k + s];
        r->assignment[pattern * n + root] = s;
      }
    if (!isfinite (optimum)) {
      status = CRPX_ERROR_OVERFLOW;
      goto done;
    }
    r->pattern_scores[pattern] = optimum;
    weighted += (long double) weights[pattern] * optimum;
    if (weighted > DBL_MAX || (mode != CRPX_PARSIMONY_CUSTOM && weighted > 9007199254740992.0L)) {
      status = CRPX_ERROR_OVERFLOW;
      goto done;
    }
    for (size_t pos = 0; pos < n; ++pos) {
      crpx_node_id node = pre[pos];
      for (size_t s = 0; s < k; ++s) {
        double candidate = crpx_private_parsimony_add (down[(size_t) node * k + s], up[(size_t) node * k + s]);
        double scale = fmax (1, fmax (fabs (candidate), fabs (optimum)));
        r->optimal[(pattern * n + node) * k + s] =
          isfinite (candidate) && fabs (candidate - optimum) <= tolerance * scale;
        total[s] = 0;
        w->zeros[s] = 0;
      }
      for (crpx_node_id child = crpx_tree_first_child (tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling (tree, child))
        for (size_t s = 0; s < k; ++s) {
          double v = message[(size_t) child * k + s];
          if (!isfinite (v))
            ++w->zeros[s];
          else
            total[s] = crpx_private_parsimony_add (total[s], v);
        }
      for (crpx_node_id child = crpx_tree_first_child (tree, node); child != CRPX_NODE_NONE;
           child = crpx_tree_next_sibling (tree, child)) {
        for (size_t s = 0; s < k; ++s) {
          double own = message[(size_t) child * k + s];
          size_t zeros = w->zeros[s] - (!isfinite (own));
          scratch[s] =
            zeros ? INFINITY
                  : crpx_private_parsimony_add (up[(size_t) node * k + s], total[s] - (isfinite (own) ? own : 0));
        }
        crpx_private_parsimony_transform (scratch, up + (size_t) child * k, k, mode, costs, true);
        size_t assigned = r->assignment[pattern * n + node], selected = 0;
        double best = INFINITY;
        for (size_t s = 0; s < k; ++s) {
          double candidate =
            crpx_private_parsimony_add (down[(size_t) child * k + s], change (assigned, s, mode, costs, k));
          if (candidate < best) {
            best = candidate;
            selected = s;
          }
        }
        r->assignment[pattern * n + child] = selected;
      }
    }
  }
  r->score = (double) weighted;
  *out = r;
  r = NULL;
done:
  crpx_parsimony_reconstruction_free (r);
  crpx_tree_workspace_reset (w->traversal);
  return status;
}

void
crpx_joint_reconstruction_free (crpx_joint_reconstruction *result)
{
  if (result) {
    free (result->assignment);
    free (result->pattern_log_joint);
    free (result);
  }
}
crpx_status
crpx_likelihood_joint_states (const crpx_tree *tree, const crpx_character_data *data,
                              const crpx_substitution_model *model, const crpx_likelihood_options *options, bool sample,
                              crpx_uniform_function uniform, void *state, crpx_reconstruction_workspace *workspace,
                              crpx_joint_reconstruction **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !data || !model || !workspace || (sample && !uniform)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t nodes = crpx_tree_node_count (tree), patterns = crpx_character_data_patterns (data);
  size_t states = crpx_substitution_model_states (model);
  if (!nodes || (patterns && nodes > SIZE_MAX / patterns / sizeof (size_t))) return CRPX_ERROR_OVERFLOW;
  crpx_status status = prepare (tree, data, workspace, true);
  crpx_joint_reconstruction *result = NULL;
  if (status != CRPX_OK) goto done;
  if (!workspace->likelihood) status = crpx_likelihood_workspace_new (&workspace->likelihood);
  if (status != CRPX_OK) goto done;
  crpx_likelihood_options settings;
  crpx_likelihood_options_init (&settings);
  if (options) settings = *options;
  settings.retain_partials = true;
  double score;
  if (sample)
    status = crpx_likelihood_evaluate (tree, data, model, &settings, workspace->likelihood, &score);
  else
    status = crpx_private_likelihood_evaluate_joint (tree, data, model, &settings, workspace->likelihood, &score);
  if (status != CRPX_OK) goto done;
  result = calloc (1, sizeof (*result));
  if (!result) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  result->nodes = nodes;
  result->patterns = patterns;
  result->assignment = patterns ? malloc (patterns * nodes * sizeof (size_t)) : NULL;
  result->pattern_log_joint = patterns ? malloc (patterns * sizeof (double)) : NULL;
  if (patterns && (!result->assignment || !result->pattern_log_joint)) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  const crpx_log_kernel *kernels = crpx_private_likelihood_kernels (workspace->likelihood);
  const double *frequencies = crpx_substitution_model_frequencies (model);
  const crpx_node_id *order = crpx_tree_workspace_preorder (workspace->traversal);
  double *weights = workspace->buffers;
  for (size_t pattern = 0; pattern < patterns; pattern++) {
    double joint = 0;
    for (size_t position = 0; position < nodes; position++) {
      crpx_node_id node = order[position], parent = crpx_tree_parent (tree, node);
      const double *partial;
      double scale;
      status = crpx_likelihood_workspace_log_partial (workspace->likelihood, node, pattern, &partial, &scale);
      if (status != CRPX_OK) goto done;
      for (size_t candidate = 0; candidate < states; candidate++)
        weights[candidate] =
          partial[candidate] + (parent == CRPX_NODE_NONE
                                  ? log (frequencies[candidate])
                                  : crpx_private_likelihood_log_transition (
                                      model, kernels + node, result->assignment[pattern * nodes + parent], candidate));
      size_t selected = 0;
      if (sample)
        status = crpx_sample_log_weights (weights, states, uniform, state, &selected);
      else {
        for (size_t candidate = 1; candidate < states; candidate++)
          if (weights[candidate] > weights[selected]) selected = candidate;
        if (!isfinite (weights[selected])) status = CRPX_ERROR_INVALID_ARGUMENT;
      }
      if (status != CRPX_OK) goto done;
      result->assignment[pattern * nodes + node] = selected;
      joint += parent == CRPX_NODE_NONE
                 ? log (frequencies[selected])
                 : crpx_private_likelihood_log_transition (model, kernels + node,
                                                           result->assignment[pattern * nodes + parent], selected);
      /* Leaf partials retain the observation weights, including soft observations. */
      if (crpx_tree_first_child (tree, node) == CRPX_NODE_NONE) joint += partial[selected] + scale;
    }
    result->pattern_log_joint[pattern] = joint;
  }
  *out = result;
  result = NULL;
done:
  crpx_joint_reconstruction_free (result);
  crpx_tree_workspace_reset (workspace->traversal);
  return status;
}
