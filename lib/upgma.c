/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * Reimplementation of biomcmc-lib distance-based tree construction.
 * BioNJ equations follow Gascuel (1997), doi:10.1093/oxfordjournals.molbev.a025808.
 */
#include "upgma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct { crpx_node_id first, second; } clustering_children;
typedef struct { crpx_node_id source, parent; } clustering_frame;

static crpx_status
initialise_distances(const struct crpx_distance_matrix_struct *distances, const char *const *labels,
                     size_t *size, double **matrix)
{
  size_t i, j, n;
  double *values, value;
  crpx_status status;
  if (!distances) return CRPX_ERROR_INVALID_ARGUMENT;
  n = crpx_distance_matrix_size(distances);
  if (!n) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n - 1 > (UINT32_MAX - 1) / 2 || n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(*values))
    return CRPX_ERROR_OVERFLOW;
  if (labels) for (i = 0; i < n; ++i) if (!labels[i]) return CRPX_ERROR_INVALID_ARGUMENT;
  values = (double *) calloc(n * n, sizeof(*values));
  if (!values) return CRPX_ERROR_ALLOCATION;
  for (i = 0; i < n; ++i) for (j = 0; j < i; ++j) {
    status = crpx_distance_matrix_get(distances, i, j, &value);
    if (status != CRPX_OK || !isfinite(value) || value < 0) {
      free(values);
      return status == CRPX_OK ? CRPX_ERROR_INVALID_ARGUMENT : status;
    }
    values[i * n + j] = values[j * n + i] = value;
  }
  *size = n;
  *matrix = values;
  return CRPX_OK;
}

static crpx_status
build_cluster_tree(size_t leaves, const char *const *labels, const clustering_children *children,
                   const double *lengths, crpx_tree **result)
{
  crpx_tree *tree = NULL;
  clustering_frame *stack, frame;
  size_t count = leaves * 2 - 1, pending;
  crpx_node_id node;
  crpx_status status;
  char label[32];
  if (count > SIZE_MAX / sizeof(*stack)) return CRPX_ERROR_OVERFLOW;
  stack = (clustering_frame *) malloc(count * sizeof(*stack));
  if (!stack) return CRPX_ERROR_ALLOCATION;
  status = crpx_tree_new(&tree);
  if (status != CRPX_OK) goto finished;
  status = crpx_tree_reserve(tree, count);
  if (status != CRPX_OK) goto finished;
  stack[0].source = (crpx_node_id) (count - 1); stack[0].parent = CRPX_NODE_NONE;
  pending = 1;
  while (pending) {
    frame = stack[--pending];
    status = crpx_tree_add_node(tree, frame.parent, &node);
    if (status != CRPX_OK) goto finished;
    if (frame.parent != CRPX_NODE_NONE) {
      status = crpx_tree_set_branch_length(tree, node, lengths[frame.source]);
      if (status != CRPX_OK) goto finished;
    }
    if (frame.source < leaves) {
      if (!labels) snprintf(label, sizeof(label), "%u", (unsigned) frame.source);
      status = crpx_tree_set_label(tree, node, labels ? labels[frame.source] : label);
      if (status != CRPX_OK) goto finished;
    }
    else {
      stack[pending].source = children[frame.source].second; stack[pending++].parent = node;
      stack[pending].source = children[frame.source].first; stack[pending++].parent = node;
    }
  }
finished:
  free(stack);
  if (status != CRPX_OK) crpx_tree_free(tree);
  else *result = tree;
  return status;
}

static crpx_status
chain_clustering(const struct crpx_distance_matrix_struct *distances, const char *const *labels,
                  bool single_linkage, crpx_tree **result)
{
  double *matrix = NULL, *height = NULL, *lengths = NULL;
  size_t *weight = NULL, *chain = NULL;
  crpx_node_id *ids = NULL;
  clustering_children *children = NULL;
  uint8_t *alive = NULL;
  size_t n, count, active, next, depth = 0, i, first, second, neighbour, previous;
  double best, value, new_height, coefficient;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  status = initialise_distances(distances, labels, &n, &matrix);
  if (status != CRPX_OK) return status;
  count = n * 2 - 1;
  if (count > SIZE_MAX / sizeof(*children) || count > SIZE_MAX / sizeof(*lengths) ||
      n > SIZE_MAX / sizeof(*chain)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
  height = (double *) calloc(n, sizeof(*height));
  lengths = (double *) calloc(count, sizeof(*lengths));
  weight = (size_t *) malloc(n * sizeof(*weight));
  chain = (size_t *) malloc(n * sizeof(*chain));
  ids = (crpx_node_id *) malloc(n * sizeof(*ids));
  children = (clustering_children *) malloc(count * sizeof(*children));
  alive = (uint8_t *) malloc(n * sizeof(*alive));
  if (!height || !lengths || !weight || !chain || !ids || !children || !alive) {
    status = CRPX_ERROR_ALLOCATION; goto finished;
  }
  for (i = 0; i < n; ++i) { weight[i] = 1; ids[i] = (crpx_node_id) i; alive[i] = 1; }
  next = n;
  active = n;
  while (active > 1) {
    if (!depth) { for (i = 0; !alive[i]; ++i) {} chain[depth++] = i; }
    first = chain[depth - 1];
    neighbour = SIZE_MAX;
    best = INFINITY;
    for (i = 0; i < n; ++i) if (alive[i] && i != first && matrix[first * n + i] < best) {
      best = matrix[first * n + i]; neighbour = i;
    }
    previous = depth > 1 ? chain[depth - 2] : SIZE_MAX;
    if (previous != SIZE_MAX && matrix[first * n + previous] <= best) neighbour = previous;
    if (neighbour != previous) { chain[depth++] = neighbour; continue; }
    second = previous;
    depth -= 2;
    new_height = matrix[first * n + second] * 0.5;
    children[next].first = ids[second]; children[next].second = ids[first];
    lengths[ids[first]] = fmax(0, new_height - height[first]);
    lengths[ids[second]] = fmax(0, new_height - height[second]);
    coefficient = (double) weight[first] / (double) (weight[first] + weight[second]);
    for (i = 0; i < n; ++i) if (alive[i] && i != first && i != second) {
      value = single_linkage ? fmin(matrix[first * n + i], matrix[second * n + i]) :
              coefficient * matrix[first * n + i] + (1 - coefficient) * matrix[second * n + i];
      matrix[second * n + i] = matrix[i * n + second] = value;
    }
    weight[second] += weight[first];
    height[second] = new_height;
    ids[second] = (crpx_node_id) next++;
    alive[first] = 0;
    --active;
  }
  status = build_cluster_tree(n, labels, children, lengths, result);
finished:
  free(matrix); free(height); free(lengths); free(weight); free(chain); free(ids); free(children); free(alive);
  return status;
}

crpx_status
crpx_tree_upgma(const struct crpx_distance_matrix_struct *distances, const char *const *labels, crpx_tree **result)
{ return chain_clustering(distances, labels, false, result); }

crpx_status
crpx_tree_single_linkage(const struct crpx_distance_matrix_struct *distances,
                         const char *const *labels, crpx_tree **result)
{ return chain_clustering(distances, labels, true, result); }

static void
store_limbs(double *lengths, crpx_node_id first, crpx_node_id second, double first_limb, double second_limb,
            double pair_distance, bool preserve_negative)
{
  if (!preserve_negative) {
    if (first_limb < 0) { first_limb = 0; second_limb = fmax(0, pair_distance); }
    else if (second_limb < 0) { second_limb = 0; first_limb = fmax(0, pair_distance); }
  }
  lengths[first] = first_limb;
  lengths[second] = second_limb;
}

static crpx_status
neighbour_joining(const struct crpx_distance_matrix_struct *distances, const char *const *labels,
                  bool bio, bool preserve_negative, crpx_tree **result)
{
  double *matrix = NULL, *variance = NULL, *lengths = NULL;
  long double *sums = NULL, criterion, best, first_limb, second_limb, lambda, delta, vij, value;
  size_t n, count, active, next, i, j, first = 0, second = 0, a, b, k;
  size_t *indices = NULL;
  crpx_node_id *ids = NULL;
  clustering_children *children = NULL;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  status = initialise_distances(distances, labels, &n, &matrix);
  if (status != CRPX_OK) return status;
  count = n * 2 - 1;
  if (count > SIZE_MAX / sizeof(*children) || count > SIZE_MAX / sizeof(*lengths) ||
      n > SIZE_MAX / sizeof(*sums) || n > SIZE_MAX / sizeof(*indices)) {
    status = CRPX_ERROR_OVERFLOW; goto finished;
  }
  if (bio) variance = (double *) malloc(n * n * sizeof(*variance));
  lengths = (double *) calloc(count, sizeof(*lengths));
  sums = (long double *) malloc(n * sizeof(*sums));
  indices = (size_t *) malloc(n * sizeof(*indices));
  ids = (crpx_node_id *) malloc(n * sizeof(*ids));
  children = (clustering_children *) malloc(count * sizeof(*children));
  if ((bio && !variance) || !lengths || !sums || !indices || !ids || !children) {
    status = CRPX_ERROR_ALLOCATION; goto finished;
  }
  for (i = 0; i < n; ++i) { indices[i] = i; ids[i] = (crpx_node_id) i; }
  if (bio) for (i = 0; i < n * n; ++i) variance[i] = matrix[i];
  active = n;
  next = n;
  while (active > 2) {
    for (i = 0; i < active; ++i) {
      sums[i] = 0;
      for (j = 0; j < active; ++j) if (j != i) sums[i] += matrix[indices[i] * n + indices[j]];
      if (!isfinite(sums[i])) { status = CRPX_ERROR_OVERFLOW; goto finished; }
    }
    best = INFINITY;
    for (i = 0; i < active; ++i) for (j = i + 1; j < active; ++j) {
      criterion = (long double) (active - 2) * matrix[indices[i] * n + indices[j]] - sums[i] - sums[j];
      if (!isfinite(criterion)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
      if (criterion < best) { best = criterion; first = i; second = j; }
    }
    a = indices[first]; b = indices[second];
    first_limb = 0.5L * ((long double) matrix[a * n + b] + (sums[first] - sums[second]) / (active - 2));
    second_limb = (long double) matrix[a * n + b] - first_limb;
    if (!isfinite((double) first_limb) || !isfinite((double) second_limb)) {
      status = CRPX_ERROR_OVERFLOW; goto finished;
    }
    lambda = 0.5L;
    vij = 0;
    if (bio) {
      vij = variance[a * n + b];
      if (vij > 0) {
        delta = 0;
        for (i = 0; i < active; ++i) if (i != first && i != second) {
          k = indices[i];
          delta += (long double) variance[b * n + k] - variance[a * n + k];
        }
        lambda += delta / (2.0L * (active - 2) * vij);
        lambda = fminl(1, fmaxl(0, lambda));
      }
    }
    children[next].first = ids[first]; children[next].second = ids[second];
    store_limbs(lengths, ids[first], ids[second], (double) first_limb, (double) second_limb,
                matrix[a * n + b], preserve_negative);
    for (i = 0; i < active; ++i) if (i != first && i != second) {
      k = indices[i];
      value = lambda * ((long double) matrix[a * n + k] - first_limb) +
              (1 - lambda) * ((long double) matrix[b * n + k] - second_limb);
      if (!isfinite((double) value)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
      matrix[a * n + k] = matrix[k * n + a] = (double) value;
      if (bio) {
        value = lambda * variance[a * n + k] + (1 - lambda) * variance[b * n + k] - lambda * (1 - lambda) * vij;
        if (!isfinite((double) value)) { status = CRPX_ERROR_OVERFLOW; goto finished; }
        variance[a * n + k] = variance[k * n + a] = fmax(0, (double) value);
      }
    }
    ids[first] = (crpx_node_id) next++;
    --active;
    indices[second] = indices[active]; ids[second] = ids[active];
  }
  if (active == 2) {
    children[next].first = ids[0]; children[next].second = ids[1];
    first_limb = (long double) matrix[indices[0] * n + indices[1]] * 0.5L;
    store_limbs(lengths, ids[0], ids[1], (double) first_limb, (double) first_limb,
                matrix[indices[0] * n + indices[1]], preserve_negative);
  }
  status = build_cluster_tree(n, labels, children, lengths, result);
finished:
  free(matrix); free(variance); free(lengths); free(sums); free(indices); free(ids); free(children);
  return status;
}

crpx_status
crpx_tree_neighbour_joining(const struct crpx_distance_matrix_struct *distances,
                            const char *const *labels, bool preserve_negative, crpx_tree **result)
{ return neighbour_joining(distances, labels, false, preserve_negative, result); }

crpx_status
crpx_tree_bionj(const struct crpx_distance_matrix_struct *distances,
                const char *const *labels, bool preserve_negative, crpx_tree **result)
{ return neighbour_joining(distances, labels, true, preserve_negative, result); }
