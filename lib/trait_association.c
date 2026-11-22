/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "trait_association.h"
#include "gaussian_internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
typedef struct {
  size_t pairs, minimum, maximum;
  bool valid;
} pairing;
struct crpx_association_workspace {
  crpx_tree_workspace *traversal;
  pairing *pairings;
  double *means, *variances;
  unsigned char *seen;
  size_t capacity;
};
static crpx_status prepare (const crpx_tree *tree, const crpx_node_id *tips, size_t count,
                            crpx_association_workspace *workspace);
static void keep (pairing *best, pairing candidate);
crpx_status
crpx_association_workspace_new (crpx_association_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void
crpx_association_workspace_free (crpx_association_workspace *workspace)
{
  if (workspace) {
    crpx_tree_workspace_free (workspace->traversal);
    free (workspace->pairings);
    free (workspace->means);
    free (workspace->variances);
    free (workspace->seen);
    free (workspace);
  }
}
crpx_status
crpx_trait_contrasting_pairs (const crpx_tree *tree, const crpx_node_id *tips, size_t count, const uint8_t *first,
                              const uint8_t *second, crpx_association_workspace *workspace, crpx_pairing_result *result)
{ /* Five states retain one exposed endpoint category or no exposed endpoint. */
  if (!first || !second || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare (tree, tips, count, workspace);
  if (status != CRPX_OK) return status;
  size_t nodes = crpx_tree_node_count (tree), observed = 0;
  memset (workspace->pairings, 0, 5 * nodes * sizeof (pairing));
  for (size_t i = 0; i < count; i++) {
    if ((first[i] > 1 && first[i] != UINT8_MAX) || (second[i] > 1 && second[i] != UINT8_MAX)) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    pairing *states = workspace->pairings + 5 * tips[i];
    states[4] = (pairing){0, 0, 0, true};
    if (first[i] != UINT8_MAX && second[i] != UINT8_MAX) {
      states[2 * first[i] + second[i]] = (pairing){0, 0, 0, true};
      observed++;
    }
  }
  const crpx_node_id *order = crpx_tree_workspace_postorder (workspace->traversal);
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = order[i], left = crpx_tree_first_child (tree, node);
    if (left == CRPX_NODE_NONE) continue;
    crpx_node_id right = crpx_tree_next_sibling (tree, left);
    pairing *target = workspace->pairings + 5 * node;
    for (size_t a = 0; a < 5; a++)
      for (size_t b = 0; b < 5; b++) {
        pairing first_pairing = workspace->pairings[5 * left + a], second_pairing = workspace->pairings[5 * right + b];
        if (!first_pairing.valid || !second_pairing.valid) continue;
        pairing candidate = {first_pairing.pairs + second_pairing.pairs, first_pairing.minimum + second_pairing.minimum,
                             first_pairing.maximum + second_pairing.maximum, true};
        if (a == 4 || b == 4)
          keep (target + (a == 4 ? b : a), candidate);
        else if ((a ^ b) == 3) {
          candidate.pairs++;
          candidate.minimum += a == 0 || a == 3;
          candidate.maximum += a == 0 || a == 3;
          keep (target + 4, candidate);
        }
      }
  }
  pairing best = {0};
  for (size_t category = 0; category < 5; category++)
    keep (&best, workspace->pairings[5 * crpx_tree_root (tree) + category]);
  *result = (crpx_pairing_result){best.pairs, best.minimum, best.maximum, observed};
done:
  crpx_tree_workspace_reset (workspace->traversal);
  return status;
}
crpx_status
crpx_trait_binary_dispersion (const crpx_tree *tree, const crpx_node_id *tips, size_t count, const uint8_t *values,
                              crpx_association_workspace *workspace, double *score)
{
  /* Propagate Gaussian nodal means, accumulating the binary-trait dispersion at each merger. */
  if (!values || !score) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = prepare (tree, tips, count, workspace);
  if (status != CRPX_OK) return status;
  size_t nodes = crpx_tree_node_count (tree);
  memset (workspace->variances, 0, nodes * sizeof (double));
  for (size_t i = 0; i < count; i++) {
    if (values[i] > 1) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    workspace->means[tips[i]] = values[i];
  }
  const crpx_node_id *order = crpx_tree_workspace_postorder (workspace->traversal);
  long double total = 0;
  for (size_t i = 0; i < nodes; i++) {
    crpx_node_id node = order[i], first = crpx_tree_first_child (tree, node);
    if (first == CRPX_NODE_NONE) continue;
    crpx_node_id second = crpx_tree_next_sibling (tree, first);
    double first_length, second_length;
    status = crpx_tree_branch_length (tree, first, &first_length);
    if (status == CRPX_OK) status = crpx_tree_branch_length (tree, second, &second_length);
    if (status != CRPX_OK) goto done;
    if (first_length <= 0 || second_length <= 0) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
    if (!isfinite (first_length + workspace->variances[first]) ||
        !isfinite (second_length + workspace->variances[second])) {
      status = CRPX_ERROR_OVERFLOW;
      goto done;
    }
    crpx_private_gaussian_product (workspace->means[first], first_length + workspace->variances[first],
                                   workspace->means[second], second_length + workspace->variances[second],
                                   workspace->means + node, workspace->variances + node);
    total += fabs (workspace->means[first] - workspace->means[node]) +
             fabs (workspace->means[second] - workspace->means[node]);
  }
  if (!isfinite (total))
    status = CRPX_ERROR_OVERFLOW;
  else
    *score = (double) total;
done:
  crpx_tree_workspace_reset (workspace->traversal);
  return status;
}
crpx_status
crpx_trait_d_statistic (double observed, double brownian_mean, double random_mean, double *value)
{
  if (!value || !isfinite (observed) || !isfinite (brownian_mean) || !isfinite (random_mean) ||
      random_mean == brownian_mean)
    return CRPX_ERROR_INVALID_ARGUMENT;
  double result = (observed - brownian_mean) / (random_mean - brownian_mean);
  if (!isfinite (result)) return CRPX_ERROR_OVERFLOW;
  *value = result;
  return CRPX_OK;
}
static void
keep (pairing *best, pairing candidate)
{
  if (!candidate.valid) return;
  if (!best->valid || candidate.pairs > best->pairs)
    *best = candidate;
  else if (candidate.pairs == best->pairs) {
    if (candidate.minimum < best->minimum) best->minimum = candidate.minimum;
    if (candidate.maximum > best->maximum) best->maximum = candidate.maximum;
  }
}
static crpx_status
prepare (const crpx_tree *tree, const crpx_node_id *tips, size_t count, crpx_association_workspace *workspace)
{
  if (!tree || !tips || !workspace || count != crpx_tree_leaf_count (tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy (tree, CRPX_TREE_REQUIRE_BINARY);
  if (status != CRPX_OK) return status;
  size_t nodes = crpx_tree_node_count (tree);
  if (!nodes) return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > SIZE_MAX / 5 / sizeof (pairing)) return CRPX_ERROR_OVERFLOW;
  if (nodes > workspace->capacity) {
    pairing *pairs = malloc (5 * nodes * sizeof (*pairs));
    double *means = malloc (nodes * sizeof (*means)), *variances = malloc (nodes * sizeof (*variances));
    unsigned char *seen = malloc (nodes);
    if (!pairs || !means || !variances || !seen) {
      free (pairs);
      free (means);
      free (variances);
      free (seen);
      return CRPX_ERROR_ALLOCATION;
    }
    free (workspace->pairings);
    free (workspace->means);
    free (workspace->variances);
    free (workspace->seen);
    workspace->pairings = pairs;
    workspace->means = means;
    workspace->variances = variances;
    workspace->seen = seen;
    workspace->capacity = nodes;
  }
  memset (workspace->seen, 0, nodes);
  for (size_t i = 0; i < count; i++) {
    if (tips[i] >= nodes || crpx_tree_first_child (tree, tips[i]) != CRPX_NODE_NONE || workspace->seen[tips[i]])
      return CRPX_ERROR_INVALID_ARGUMENT;
    workspace->seen[tips[i]] = 1;
  }
  if (!workspace->traversal) status = crpx_tree_workspace_new (CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace->traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh (workspace->traversal, tree);
  return status;
}
