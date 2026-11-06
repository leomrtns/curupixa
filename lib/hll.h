/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_HLL_H
#define CURUPIXA_HLL_H
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_hll crpx_hll;
/* Classic 64-bit HyperLogLog; precision 4..18. A mutable sketch belongs to one worker;
 * merge worker/rank sketches after synchronisation. Seeds and precision must match. */
crpx_status crpx_hll_new(unsigned precision, uint64_t seed, crpx_hll **result);
void crpx_hll_free(crpx_hll *sketch);
void crpx_hll_reset(crpx_hll *sketch);
crpx_status crpx_hll_clone(const crpx_hll *source, crpx_hll **result);
crpx_status crpx_hll_add(crpx_hll *sketch, const void *data, size_t length);
/* Prehashed values must have uniformly distributed bits, with one common hash convention. */
crpx_status crpx_hll_add_hash(crpx_hll *sketch, uint64_t hash);
crpx_status crpx_hll_merge(crpx_hll *destination, const crpx_hll *source);
crpx_status crpx_hll_estimate(const crpx_hll *sketch, double *estimate, double *relative_standard_error);
size_t crpx_hll_register_count(const crpx_hll *sketch);
/* Stable little-endian wire format, suitable for files or MPI byte messages. Calling with
 * buffer NULL reports the required size through written and returns OUT_OF_RANGE. */
crpx_status crpx_hll_serialize(const crpx_hll *sketch, void *buffer, size_t capacity, size_t *written);
crpx_status crpx_hll_deserialize(const void *buffer, size_t length, crpx_hll **result);
#ifdef __cplusplus
}
#endif
#endif
