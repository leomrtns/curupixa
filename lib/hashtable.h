/* Reimplementation of biomcmc-lib hash containers (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_HASHTABLE_H
#define CURUPIXA_HASHTABLE_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_hashtable_struct *crpx_hashtable_t;
/* Byte-key -> uint64_t map; keys are copied, including embedded NULs. Seed is
 * explicit and deterministic. Concurrent reads are safe between mutations.
 * This non-cryptographic hash is intended for scientific data. */
crpx_status crpx_hashtable_new(size_t expected_items, uint64_t seed, crpx_hashtable_t *out);
void crpx_hashtable_free(crpx_hashtable_t table);
void crpx_hashtable_clear(crpx_hashtable_t table);
size_t crpx_hashtable_size(const struct crpx_hashtable_struct *table);
crpx_status crpx_hashtable_put(crpx_hashtable_t table, const void *key, size_t length, uint64_t value,
                             bool *replaced);
crpx_status crpx_hashtable_get(const struct crpx_hashtable_struct *table, const void *key, size_t length,
                             uint64_t *value);
crpx_status crpx_hashtable_remove(crpx_hashtable_t table, const void *key, size_t length);
/* cursor must initially be zero. Key remains owned by the table and valid until
 * mutation/free; iteration order is unspecified. NOT_FOUND means end. */
crpx_status crpx_hashtable_next(const struct crpx_hashtable_struct *table, size_t *cursor,
                              const void **key, size_t *length, uint64_t *value);
crpx_status crpx_hashtable_put_string(crpx_hashtable_t table, const char *key, uint64_t value, bool *replaced);
crpx_status crpx_hashtable_get_string(const struct crpx_hashtable_struct *table, const char *key, uint64_t *value);

#ifdef __cplusplus
}
#endif
#endif
