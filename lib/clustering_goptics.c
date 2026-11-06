/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "clustering_goptics.h"
#include "quickselect_quantile.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

struct crpx_optics_struct {
  size_t size;
  double epsilon;
  size_t *order, *predecessors;
  double *core, *reachability;
  size_t *border_predecessors;
  double *border_reachability;
};

void crpx_optics_free(crpx_optics_t ordering)
{
  if (ordering) {
    free(ordering->order); free(ordering->predecessors);
    free(ordering->core); free(ordering->reachability);
    free(ordering->border_predecessors); free(ordering->border_reachability); free(ordering);
  }
}
size_t crpx_optics_size(const struct crpx_optics_struct *ordering) { return ordering ? ordering->size : 0; }
const size_t *crpx_optics_order(const struct crpx_optics_struct *ordering) { return ordering ? ordering->order : NULL; }
const double *crpx_optics_core_distances(const struct crpx_optics_struct *ordering) { return ordering ? ordering->core : NULL; }
const double *crpx_optics_reachability_distances(const struct crpx_optics_struct *ordering)
{ return ordering ? ordering->reachability : NULL; }
const size_t *crpx_optics_predecessors(const struct crpx_optics_struct *ordering)
{ return ordering ? ordering->predecessors : NULL; }

crpx_status crpx_optics_run(crpx_distance_generator_t generator, size_t component, size_t min_points,
                           double epsilon, crpx_optics_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t size = crpx_distance_generator_samples(generator);
  if (!generator || !size || !min_points || min_points > size || component >= crpx_distance_generator_components(generator) ||
      !isfinite(epsilon) || epsilon < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof(double) || size > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  crpx_optics_t ordering = calloc(1, sizeof(*ordering));
  if (!ordering) return CRPX_ERROR_ALLOCATION;
  ordering->size = size; ordering->epsilon = epsilon;
  ordering->order = malloc(size * sizeof(size_t));
  ordering->predecessors = malloc(size * sizeof(size_t));
  ordering->core = malloc(size * sizeof(double));
  ordering->reachability = malloc(size * sizeof(double));
  ordering->border_predecessors = malloc(size * sizeof(size_t));
  ordering->border_reachability = malloc(size * sizeof(double));
  bool *processed = calloc(size, sizeof(bool));
  double *distances = malloc(size * sizeof(double)), *neighbours = malloc(size * sizeof(double));
  if (!ordering->order || !ordering->predecessors || !ordering->core || !ordering->reachability ||
      !ordering->border_predecessors || !ordering->border_reachability || !processed || !distances || !neighbours) {
    free(processed); free(distances); free(neighbours); crpx_optics_free(ordering); return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < size; ++i) {
    ordering->core[i] = ordering->reachability[i] = INFINITY;
    ordering->predecessors[i] = SIZE_MAX;
    ordering->border_predecessors[i] = SIZE_MAX;
    ordering->border_reachability[i] = INFINITY;
  }
  crpx_status status = CRPX_OK;
  for (size_t position = 0; position < size; ++position) {
    size_t current = SIZE_MAX;
    for (size_t i = 0; i < size; ++i)
      if (!processed[i] && (current == SIZE_MAX || ordering->reachability[i] < ordering->reachability[current])) current = i;
    processed[current] = true;
    ordering->order[position] = current;
    size_t count = 0;
    for (size_t i = 0; i < size; ++i) {
      status = crpx_distance_generator_get(generator, current, i, component, distances + i);
      if (status != CRPX_OK) goto done;
      if (distances[i] < 0 || !isfinite(distances[i])) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      if (distances[i] <= epsilon) neighbours[count++] = distances[i];
    }
    if (count >= min_points) {
      status = crpx_select_double(neighbours, count, min_points - 1, ordering->core + current);
      if (status != CRPX_OK) goto done;
      for (size_t i = 0; i < size; ++i)
        if (distances[i] <= epsilon) {
          double reachability = fmax(ordering->core[current], distances[i]);
          if (reachability < ordering->border_reachability[i] ||
              (reachability == ordering->border_reachability[i] && current < ordering->border_predecessors[i])) {
            ordering->border_reachability[i] = reachability;
            ordering->border_predecessors[i] = current;
          }
          if (!processed[i] && reachability < ordering->reachability[i]) {
            ordering->reachability[i] = reachability;
            ordering->predecessors[i] = current;
          }
        }
    }
  }
done:
  free(processed); free(distances); free(neighbours);
  if (status != CRPX_OK) { crpx_optics_free(ordering); return status; }
  *out = ordering;
  return CRPX_OK;
}

crpx_status crpx_optics_extract_clusters(const struct crpx_optics_struct *ordering, double epsilon,
                                        size_t *labels, size_t capacity, size_t *cluster_count)
{
  if (!ordering || !labels || !cluster_count || !isfinite(epsilon) || epsilon < 0 || epsilon > ordering->epsilon)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity < ordering->size) return CRPX_ERROR_OUT_OF_RANGE;
  size_t count = 0, cluster = SIZE_MAX;
  for (size_t position = 0; position < ordering->size; ++position) {
    size_t sample = ordering->order[position];
    if (ordering->reachability[sample] > epsilon) {
      if (ordering->core[sample] <= epsilon) { cluster = count++; labels[sample] = cluster; }
      else labels[sample] = SIZE_MAX;
    }
    else labels[sample] = cluster;
  }
  /* A border point may precede the first core in its component and therefore
   * retain infinite OPTICS reachability. Attach such points after core cluster
   * labels are known, without altering the ordering's reachability semantics. */
  for (size_t sample = 0; sample < ordering->size; ++sample)
    if (labels[sample] == SIZE_MAX && ordering->border_reachability[sample] <= epsilon) {
      size_t core = ordering->border_predecessors[sample];
      if (core != SIZE_MAX) labels[sample] = labels[core];
    }
  *cluster_count = count;
  return CRPX_OK;
}
