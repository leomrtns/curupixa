/* Inspired by biomcmc-lib distance_generator.h.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_DISTANCE_GENERATOR_H
#define CURUPIXA_DISTANCE_GENERATOR_H
#include <stdbool.h>
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_distance_generator_struct *crpx_distance_generator_t;
typedef crpx_status (*crpx_distance_callback)(void *data, size_t first, size_t second,
                                             double *distances, size_t components);

/* The callback receives first<second and writes components finite distances.
 * Negative distances are permitted. Diagonal distances are zero without a call.
 * data remains caller-owned and must outlive use of the generator.
 * Calls/reset/replacement are serialised, including callback execution. The
 * callback must not recursively access this generator. Destruction requires
 * exclusive ownership. Cached vectors are allocated only on first use; the
 * optional pointer index requires O(samples^2) storage. cache=false uses O(1)
 * persistent storage. Different generators may share data only if it is safe. */
crpx_status crpx_distance_generator_new(size_t samples, size_t components, crpx_distance_callback callback,
                                       void *data, bool cache, crpx_distance_generator_t *out);
void crpx_distance_generator_free(crpx_distance_generator_t generator);
size_t crpx_distance_generator_samples(const struct crpx_distance_generator_struct *generator);
size_t crpx_distance_generator_components(const struct crpx_distance_generator_struct *generator);
crpx_status crpx_distance_generator_get(crpx_distance_generator_t generator, size_t first, size_t second,
                                       size_t component, double *out);
crpx_status crpx_distance_generator_get_all(crpx_distance_generator_t generator, size_t first, size_t second,
                                           double *out, size_t capacity);
crpx_status crpx_distance_generator_reset(crpx_distance_generator_t generator);
crpx_status crpx_distance_generator_set_callback(crpx_distance_generator_t generator,
                                                crpx_distance_callback callback, void *data);

#ifdef __cplusplus
}
#endif
#endif
