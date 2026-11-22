/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Cluster objects from supplied distances using the generalised OPTICS calculation. */
#include "clustering_goptics.h"
#include "quickselect_quantile.h"
#include "priority_queue_internal.h"
#include "clustering_internal.h"
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

struct optics_distances { crpx_distance_generator_t generator; size_t component; };
static crpx_status optics_neighbours (void *data, size_t query, double epsilon, size_t *ids,
                                      double *distances, size_t capacity, size_t *count);

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

crpx_status
crpx_optics_run (crpx_distance_generator_t generator, size_t component, size_t min_points,
                 double epsilon, crpx_optics_t *out)
{ // Ankerst et al., SIGMOD (1999): use the same ordering kernel for a full scan or an indexed neighbourhood search.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!generator || component >= crpx_distance_generator_components (generator)) return CRPX_ERROR_INVALID_ARGUMENT;
  struct optics_distances data = {generator, component};
  return crpx_optics_run_neighbours (crpx_distance_generator_samples (generator), optics_neighbours, &data,
                                     min_points, epsilon, out);
}

crpx_status
crpx_optics_run_neighbours (size_t size, crpx_neighbour_callback neighbours, void *data,
                            size_t min_points, double epsilon, crpx_optics_t *out)
{ // Decrease-key updates replace the previous scan of all remaining points at each ordering step.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!neighbours || !size || !min_points || min_points > size || !isfinite (epsilon) || epsilon < 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_priority_queue queue;
  crpx_status status = crpx_queue_init (&queue, size);
  if (status != CRPX_OK) return status;
  crpx_optics_t ordering = calloc (1, sizeof (*ordering));
  bool *processed = calloc (size, sizeof (*processed));
  size_t *ids = malloc (size * sizeof (*ids)), *seen = calloc (size, sizeof (*seen));
  double *distances = malloc (size * sizeof (*distances)), *selection = malloc (size * sizeof (*selection));
  if (!ordering || !processed || !ids || !seen || !distances || !selection) { status = CRPX_ERROR_ALLOCATION; goto done; }
  ordering->size = size; ordering->epsilon = epsilon;
  ordering->order = malloc (size * sizeof (*ordering->order));
  ordering->predecessors = malloc (size * sizeof (*ordering->predecessors));
  ordering->core = malloc (size * sizeof (*ordering->core));
  ordering->reachability = malloc (size * sizeof (*ordering->reachability));
  ordering->border_predecessors = malloc (size * sizeof (*ordering->border_predecessors));
  ordering->border_reachability = malloc (size * sizeof (*ordering->border_reachability));
  if (!ordering->order || !ordering->predecessors || !ordering->core || !ordering->reachability ||
      !ordering->border_predecessors || !ordering->border_reachability) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i = 0; i < size; i++) {
    ordering->core[i] = ordering->reachability[i] = ordering->border_reachability[i] = INFINITY;
    ordering->predecessors[i] = ordering->border_predecessors[i] = SIZE_MAX;
  }
  size_t root = 0;
  for (size_t position = 0; position < size; position++) {
    size_t current;
    if (queue.size) current = crpx_queue_pop (&queue);
    else { while (processed[root]) root++; current = root; }
    processed[current] = true; ordering->order[position] = current;
    size_t count = 0;
    status = neighbours (data, current, epsilon, ids, distances, size - 1, &count);
    if (status != CRPX_OK) goto done;
    if (count > size - 1) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    for (size_t i = 0; i < count; i++) {
      if (ids[i] >= size || ids[i] == current || seen[ids[i]] == position + 1 ||
          !isfinite (distances[i]) || distances[i] < 0 || distances[i] > epsilon) {
        status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
      }
      seen[ids[i]] = position + 1; selection[i] = distances[i];
    }
    if (min_points == 1) ordering->core[current] = 0;
    else if (count >= min_points - 1) {
      status = crpx_select_double (selection, count, min_points - 2, ordering->core + current);
      if (status != CRPX_OK) goto done;
    }
    if (!isfinite (ordering->core[current])) continue;
    for (size_t i = 0; i < count; i++) {
      size_t next = ids[i];
      double reachability = fmax (ordering->core[current], distances[i]);
      if (reachability < ordering->border_reachability[next] ||
          (reachability == ordering->border_reachability[next] && current < ordering->border_predecessors[next])) {
        ordering->border_reachability[next] = reachability; ordering->border_predecessors[next] = current;
      }
      if (!processed[next] && reachability < ordering->reachability[next]) {
        ordering->reachability[next] = reachability; ordering->predecessors[next] = current;
        crpx_queue_offer (&queue, next, reachability);
      }
    }
  }
done:
  free (processed); free (ids); free (seen); free (distances); free (selection); crpx_queue_free (&queue);
  if (status != CRPX_OK) crpx_optics_free (ordering);
  else *out = ordering;
  return status;
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
  /* A border point may precede the first core in its component and therefore retain infinite OPTICS reachability. Attach
   such points after core cluster labels are known, without altering the ordering's reachability semantics. */
  for (size_t sample = 0; sample < ordering->size; ++sample)
    if (labels[sample] == SIZE_MAX && ordering->border_reachability[sample] <= epsilon) {
      size_t core = ordering->border_predecessors[sample];
      if (core != SIZE_MAX) labels[sample] = labels[core];
    }
  *cluster_count = count;
  return CRPX_OK;
}

static crpx_status
optics_neighbours (void *data, size_t query, double epsilon, size_t *ids,
                    double *distances, size_t capacity, size_t *count)
{ // The default needs no geometric assumptions; an indexed callback can replace this quadratic scan.
  struct optics_distances *source = data;
  size_t size = crpx_distance_generator_samples (source->generator);
  *count = 0;
  for (size_t i = 0; i < size; i++) if (i != query) {
    double distance;
    crpx_status status = crpx_cluster_distance (source->generator, query, i, source->component, &distance);
    if (status != CRPX_OK) return status;
    if (distance <= epsilon) {
      if (*count == capacity) return CRPX_ERROR_OUT_OF_RANGE;
      ids[*count] = i; distances[(*count)++] = distance;
    }
  }
  return CRPX_OK;
}
