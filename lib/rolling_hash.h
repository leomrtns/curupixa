/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_ROLLING_HASH_H
#define CURUPIXA_ROLLING_HASH_H
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_rolling_hash crpx_rolling_hash;
/* Seeded polynomial 64-bit byte-window hash; arbitrary positive window lengths. The initial
 * window takes O(window) time and successive windows O(1). Hashes are non-cryptographic. */
crpx_status crpx_rolling_hash_new(size_t window, uint64_t seed, crpx_rolling_hash **result);
void crpx_rolling_hash_free(crpx_rolling_hash *iterator);
/* Borrowed input must remain unchanged until reset/free. Embedded zero bytes are supported. */
crpx_status crpx_rolling_hash_reset(crpx_rolling_hash *iterator, const void *data, size_t length);
crpx_status crpx_rolling_hash_next(crpx_rolling_hash *iterator, uint64_t *hash, size_t *position);
#ifdef __cplusplus
}
#endif
#endif
