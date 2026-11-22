/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "tree_prior.h"
#include "tree_comparison.h"
#include <math.h>
#include <stdlib.h>

static crpx_status validate_prior (const crpx_tree_prior *prior);
static double bd_cdf (double age, double birth, double death);
static double bd_log_density (double age, double birth, double death);
static int compare_double (const void *a, const void *b);
static crpx_status draw_index (size_t count, crpx_uniform_function uniform, void *state, size_t *out);

crpx_status
crpx_tree_prior_sample (const crpx_tree_prior *prior, size_t leaves, const char *const *labels,
                        crpx_uniform_function uniform, void *state, crpx_tree **out)
{ /* Construct a topology independently of biological waiting-time parameters. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_status status = validate_prior (prior);
  if (status != CRPX_OK) return status;
  if (!leaves || !uniform) return CRPX_ERROR_INVALID_ARGUMENT;
  if (leaves > UINT32_MAX / 2 || leaves > SIZE_MAX / 4 / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  size_t nodes = 2 * leaves - 1, root = 0;
  size_t *children = calloc (2 * leaves, sizeof (*children));
  size_t *active = malloc (nodes * sizeof (*active)), *parent = malloc (nodes * sizeof (*parent));
  double *ages = calloc (nodes, sizeof (*ages)), *times = calloc (leaves, sizeof (*times));
  if (!children || !active || !parent || !ages || !times) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  if (prior->kind == CRPX_TREE_PRIOR_UNIFORM && leaves > 1) {
    root = leaves;
    children[0] = 0;
    children[1] = 1;
    active[0] = 0;
    active[1] = 1;
    active[2] = root;
    parent[0] = parent[1] = root;
    parent[root] = SIZE_MAX;
    size_t used = 3;
    for (size_t tip = 2; tip < leaves; tip++) {
      size_t selected;
      status = draw_index (used, uniform, state, &selected);
      if (status != CRPX_OK) goto done;
      size_t edge = active[selected], old_parent = parent[edge], inserted = leaves + tip - 1;
      children[2 * (inserted - leaves)] = edge;
      children[2 * (inserted - leaves) + 1] = tip;
      parent[inserted] = old_parent;
      parent[edge] = parent[tip] = inserted;
      if (old_parent == SIZE_MAX)
        root = inserted;
      else
        children[2 * (old_parent - leaves) + (children[2 * (old_parent - leaves)] != edge)] = inserted;
      active[used++] = inserted;
      active[used++] = tip;
    }
  } else {
    double age = 0;
    double birth = prior->birth_rate, death = prior->kind == CRPX_TREE_PRIOR_YULE ? 0 : prior->death_rate;
    for (size_t i = 0; i + 1 < leaves; i++) {
      double draw;
      status = crpx_sample_uniform (uniform, state, &draw);
      if (status != CRPX_OK) goto done;
      /* A midpoint convention avoids infinite waiting times for a valid zero draw. */
      if (draw == 0) draw = 0x1p-54;
      if (prior->kind == CRPX_TREE_PRIOR_COALESCENT) {
        double lineages = (double) (leaves - i);
        age += -log1p (-draw) * prior->coalescent_scale / (lineages * (lineages - 1) / 2);
        times[i] = age;
      } else {
        double probability = draw * bd_cdf (prior->origin_age, birth, death), difference = birth - death;
        times[i] = difference == 0 ? probability / (birth * (1 - probability))
                                   : -log1p (-(difference * probability) / (birth - death * probability)) / difference;
      }
      if (!isfinite (times[i]) || times[i] <= 0) {
        status = CRPX_ERROR_OVERFLOW;
        goto done;
      }
    }
    if (prior->kind != CRPX_TREE_PRIOR_COALESCENT) qsort (times, leaves - 1, sizeof (*times), compare_double);
    for (size_t i = 0; i < leaves; i++)
      active[i] = i;
    for (size_t remaining = leaves; remaining > 1; remaining--) {
      size_t first, second, node = leaves + (leaves - remaining);
      status = draw_index (remaining, uniform, state, &first);
      if (status != CRPX_OK) goto done;
      size_t left = active[first];
      active[first] = active[remaining - 1];
      status = draw_index (remaining - 1, uniform, state, &second);
      if (status != CRPX_OK) goto done;
      children[2 * (node - leaves)] = left;
      children[2 * (node - leaves) + 1] = active[second];
      active[second] = node;
      ages[node] = times[leaves - remaining];
      root = node;
    }
  }
  status = crpx_tree_from_binary_history (leaves, children, root, prior->kind == CRPX_TREE_PRIOR_UNIFORM ? NULL : ages,
                                          labels, out);
done:
  free (children);
  free (active);
  free (parent);
  free (ages);
  free (times);
  return status;
}

crpx_status
crpx_tree_prior_log_density (const crpx_tree_prior *prior, const crpx_tree *tree, double *value)
{
  crpx_status status = validate_prior (prior);
  if (status != CRPX_OK || !tree || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  size_t leaves = crpx_tree_leaf_count (tree), nodes = crpx_tree_node_count (tree);
  if (!leaves) return CRPX_ERROR_INVALID_ARGUMENT;
  double score = 0;
  if (prior->kind == CRPX_TREE_PRIOR_UNIFORM) {
    for (size_t i = 2; i <= leaves; i++)
      score -= log (2 * (double) i - 3);
    *value = score;
    return CRPX_OK;
  }
  crpx_tree_workspace *traversal = NULL;
  double *ages = calloc (nodes, sizeof (*ages)), *times = malloc (leaves * sizeof (*times));
  if (!ages || !times) {
    free (ages);
    free (times);
    return CRPX_ERROR_ALLOCATION;
  }
  status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (traversal, tree);
  if (status != CRPX_OK) goto done;
  const crpx_node_id *order = crpx_tree_workspace_postorder (traversal);
  size_t count = 0;
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = order[i], child = crpx_tree_first_child (tree, node);
    if (child == CRPX_NODE_NONE) continue;
    bool first = true;
    for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling (tree, child)) {
      double length;
      status = crpx_tree_branch_length (tree, child, &length);
      if (status != CRPX_OK) goto done;
      if (length <= 0 || !isfinite (length)) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      double candidate = ages[child] + length;
      if (!first && fabs (candidate - ages[node]) > 1e-9 * fmax (1, candidate)) {
        status = CRPX_ERROR_INVALID_ARGUMENT;
        goto done;
      }
      ages[node] = candidate;
      first = false;
    }
    times[count++] = ages[node];
  }
  qsort (times, count, sizeof (*times), compare_double);
  if (prior->kind == CRPX_TREE_PRIOR_COALESCENT) {
    double previous = 0;
    for (size_t i = 0; i < count; i++) {
      double remaining = (double) (leaves - i), rate = remaining * (remaining - 1) / (2 * prior->coalescent_scale);
      score -= log (prior->coalescent_scale) + rate * (times[i] - previous);
      previous = times[i];
    }
  } else {
    double birth = prior->birth_rate, death = prior->kind == CRPX_TREE_PRIOR_YULE ? 0 : prior->death_rate;
    double normaliser = log (bd_cdf (prior->origin_age, birth, death));
    score = lgamma ((double) leaves);
    for (size_t i = 0; i < count; i++) {
      if (times[i] >= prior->origin_age) {
        score = -INFINITY;
        break;
      }
      double remaining = (double) (leaves - i);
      score += bd_log_density (times[i], birth, death) - normaliser - log (remaining * (remaining - 1) / 2);
    }
  }
  if (isnan (score))
    status = CRPX_ERROR_OVERFLOW;
  else
    *value = score;
done:
  free (ages);
  free (times);
  crpx_tree_workspace_free (traversal);
  return status;
}

crpx_status
crpx_tree_resolution_sample (const crpx_tree *constraint, crpx_uniform_function uniform, void *state, crpx_tree **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!constraint) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t node = 0; node < crpx_tree_node_count (constraint); node++) {
    crpx_node_id child = crpx_tree_first_child (constraint, (crpx_node_id) node);
    if (child != CRPX_NODE_NONE && crpx_tree_next_sibling (constraint, child) == CRPX_NODE_NONE)
      return CRPX_ERROR_UNSUPPORTED;
  }
  return crpx_tree_resolve_random (constraint, uniform, state, out);
}

crpx_status
crpx_tree_resolution_log_density (const crpx_tree *constraint, const crpx_node_id *constraint_tips,
                                  const crpx_tree *candidate, const crpx_node_id *candidate_tips, size_t taxa,
                                  double *value)
{
  if (!constraint || !candidate || !constraint_tips || !candidate_tips || !value ||
      taxa != crpx_tree_leaf_count (constraint) || taxa != crpx_tree_leaf_count (candidate))
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (candidate, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  crpx_tree_comparison *first = NULL, *second = NULL;
  status = crpx_tree_comparison_new (constraint, constraint_tips, taxa, true, false, &first);
  if (status == CRPX_OK) status = crpx_tree_comparison_new (candidate, candidate_tips, taxa, true, false, &second);
  crpx_tree_distance_result comparison;
  if (status == CRPX_OK) status = crpx_tree_comparison_distance (first, second, &comparison);
  crpx_tree_comparison_free (first);
  crpx_tree_comparison_free (second);
  if (status != CRPX_OK) return status;
  /* RF=a+b-2*common: every constraint clade must occur in the candidate. */
  if (comparison.second_splits < comparison.first_splits ||
      comparison.rf != comparison.second_splits - comparison.first_splits) {
    *value = -INFINITY;
    return CRPX_OK;
  }
  double score = 0;
  for (size_t node = 0; node < crpx_tree_node_count (constraint); node++) {
    size_t children = 0;
    for (crpx_node_id child = crpx_tree_first_child (constraint, (crpx_node_id) node); child != CRPX_NODE_NONE;
         child = crpx_tree_next_sibling (constraint, child))
      children++;
    if (children == 1) return CRPX_ERROR_UNSUPPORTED;
    for (size_t i = 2; i <= children; i++)
      score -= log (2 * (double) i - 3);
  }
  *value = score;
  return CRPX_OK;
}

static crpx_status
validate_prior (const crpx_tree_prior *prior)
{
  if (!prior || prior->kind < CRPX_TREE_PRIOR_UNIFORM || prior->kind > CRPX_TREE_PRIOR_COALESCENT)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (prior->kind == CRPX_TREE_PRIOR_UNIFORM) return CRPX_OK;
  if (prior->kind == CRPX_TREE_PRIOR_COALESCENT)
    return isfinite (prior->coalescent_scale) && prior->coalescent_scale > 0 ? CRPX_OK : CRPX_ERROR_INVALID_ARGUMENT;
  if (!isfinite (prior->birth_rate) || prior->birth_rate <= 0 || !isfinite (prior->origin_age) ||
      prior->origin_age <= 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (prior->kind == CRPX_TREE_PRIOR_BIRTH_DEATH && (!isfinite (prior->death_rate) || prior->death_rate < 0))
    return CRPX_ERROR_INVALID_ARGUMENT;
  return CRPX_OK;
}
static double
bd_cdf (double age, double birth, double death)
{
  double difference = birth - death;
  if (!difference) {
    double product = birth * age;
    return isinf (product) ? 1 : product / (1 + product);
  }
  if (difference < 0) {
    double loss = -expm1 (difference * age);
    return birth * loss / (-difference + birth * loss);
  }
  double z = -expm1 (-difference * age);
  return birth * z / (difference + death * z);
}
static double
bd_log_density (double age, double birth, double death)
{
  double difference = birth - death;
  if (!difference) return log (birth) - 2 * log1p (birth * age);
  if (difference < 0) {
    double denominator = -difference + birth * (-expm1 (difference * age));
    return log (birth) + 2 * log (-difference) + difference * age - 2 * log (denominator);
  }
  double denominator = difference + death * (-expm1 (-difference * age));
  return log (birth) + 2 * log (difference) - difference * age - 2 * log (denominator);
}
static int
compare_double (const void *a, const void *b)
{
  double first = *(const double *) a, second = *(const double *) b;
  return (first > second) - (first < second);
}
static crpx_status
draw_index (size_t count, crpx_uniform_function uniform, void *state, size_t *out)
{
  double draw;
  crpx_status status = crpx_sample_uniform (uniform, state, &draw);
  if (status == CRPX_OK) *out = (size_t) ((long double) draw * count);
  return status;
}
