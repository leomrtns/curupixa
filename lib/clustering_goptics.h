/* Inspired by biomcmc-lib's GOPTICS interface.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_CLUSTERING_GOPTICS_H
#define CURUPIXA_CLUSTERING_GOPTICS_H
#include <stddef.h>
#include "distance_generator.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_optics_struct *crpx_optics_t;
/* OPTICS ordering with O(n) workspace and O(n^2) distance queries. Distances must
 * be finite and non-negative. min_points includes the query point itself and
 * must lie in 1..n. Neighbourhoods include distance == epsilon. Ties use sample ID.
 * The result owns all arrays; it does not retain the generator. */
crpx_status crpx_optics_run(crpx_distance_generator_t generator, size_t component, size_t min_points,
                           double epsilon, crpx_optics_t *out);
void crpx_optics_free(crpx_optics_t ordering);
size_t crpx_optics_size(const struct crpx_optics_struct *ordering);
/* order[position] identifies an original sample. The remaining arrays use
 * original sample indices, with INFINITY for undefined distances and SIZE_MAX
 * for undefined predecessors. Concurrent result reads are safe. */
const size_t *crpx_optics_order(const struct crpx_optics_struct *ordering);
const double *crpx_optics_core_distances(const struct crpx_optics_struct *ordering);
const double *crpx_optics_reachability_distances(const struct crpx_optics_struct *ordering);
const size_t *crpx_optics_predecessors(const struct crpx_optics_struct *ordering);
/* Extract DBSCAN clusters at epsilon <= the ordering's search epsilon. labels
 * has n entries in original sample order, cluster IDs start at zero, and noise
 * is SIZE_MAX. Leading border points attach to their best adjacent core after
 * the ordering cut; ties use core sample ID. Extraction does not mutate the ordering. */
crpx_status crpx_optics_extract_clusters(const struct crpx_optics_struct *ordering, double epsilon,
                                        size_t *labels, size_t capacity, size_t *cluster_count);

#ifdef __cplusplus
}
#endif
#endif
