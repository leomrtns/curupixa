/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file kmedoids.h
    PAM medoids for arbitrary finite non-negative dissimilarities, using FastPAM1 swap evaluation.
    The result records representatives, assignments and the sum of distances to the representatives. */
#ifndef CURUPIXA_KMEDOIDS_H
#define CURUPIXA_KMEDOIDS_H
#include <stdbool.h>
#include <stddef.h>
#include "distance_generator.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  size_t samples, clusters, swaps;
  size_t *medoids, *labels;
  double *distances;
  double cost;
  bool converged;
} crpx_kmedoids;

/* initial contains k distinct sample IDs, or NULL for deterministic PAM BUILD initialisation.
   max_swaps bounds accepted swaps; zero evaluates the initial representatives. The result remains usable when the
   limit is reached (converged=false). Each FastPAM1 pass uses O(n^2+nk) queries/operations and O(n+k) scratch;
   BUILD costs O(kn^2). Equal costs use sample ID. labels index the medoids array; a medoid labels itself.
   Schubert and Rousseeuw (2021), doi:10.1016/j.is.2021.101804. */
crpx_status crpx_kmedoids_run (crpx_distance_generator_t generator, size_t component, size_t k,
                              const size_t *initial, size_t max_swaps, crpx_kmedoids **out);
void crpx_kmedoids_free (crpx_kmedoids *result);
#ifdef __cplusplus
}
#endif
#endif
