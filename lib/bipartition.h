/* Inspired by biomcmc-lib bipartition.h (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_BIPARTITION_H
#define CURUPIXA_BIPARTITION_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_bipartition_struct *crpx_bipartition_t;
typedef enum {
  CRPX_BITS_OR, CRPX_BITS_AND, CRPX_BITS_ANDNOT, CRPX_BITS_XOR, CRPX_BITS_XORNOT
} crpx_bit_operation;

/* Independent ownership; concurrent reads are safe. Serialise changes to the
 * same object. Padding bits always remain zero. Output may alias either input. */
crpx_status crpx_bipartition_new(size_t bits, crpx_bipartition_t *out);
crpx_status crpx_bipartition_clone(const struct crpx_bipartition_struct *source, crpx_bipartition_t *out);
void crpx_bipartition_free(crpx_bipartition_t bits);
size_t crpx_bipartition_size(const struct crpx_bipartition_struct *bits);
size_t crpx_bipartition_count(const struct crpx_bipartition_struct *bits);
crpx_status crpx_bipartition_resize(crpx_bipartition_t bits, size_t size);
crpx_status crpx_bipartition_set(crpx_bipartition_t bits, size_t position, bool value);
crpx_status crpx_bipartition_get(const struct crpx_bipartition_struct *bits, size_t position, bool *value);
void crpx_bipartition_clear(crpx_bipartition_t bits);
crpx_status crpx_bipartition_singleton(crpx_bipartition_t bits, size_t position);
crpx_status crpx_bipartition_copy(crpx_bipartition_t to, const struct crpx_bipartition_struct *from);
crpx_status crpx_bipartition_combine(crpx_bipartition_t out, const struct crpx_bipartition_struct *a,
                                   const struct crpx_bipartition_struct *b, crpx_bit_operation operation);
crpx_status crpx_bipartition_complement(crpx_bipartition_t out, const struct crpx_bipartition_struct *source);
bool crpx_bipartition_equal(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b);
bool crpx_bipartition_same_split(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b);
bool crpx_bipartition_contains(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b);
/* Unrooted splits are compatible if one of their four intersections is empty. */
bool crpx_bipartition_compatible(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b);
/* Choose the smaller side; break equal-size ties by clearing bit zero. */
void crpx_bipartition_canonicalise(crpx_bipartition_t bits);
crpx_status crpx_bipartition_indices(const struct crpx_bipartition_struct *bits, size_t *indices,
                                   size_t capacity, size_t *required);
/* Read-only words are invalidated by resize/free. */
const uint64_t *crpx_bipartition_words(const struct crpx_bipartition_struct *bits, size_t *count);

#ifdef __cplusplus
}
#endif
#endif
