/* Reimplementation of biomcmc-lib char_vector (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_CHAR_VECTOR_H
#define CURUPIXA_CHAR_VECTOR_H
#include <stdbool.h>
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_char_vector_struct *crpx_char_vector_t;

/* Owns copies of all strings. Empty slots read as "". No reference counting or
 * borrowed writable strings. Concurrent reads are safe until the next change. */
crpx_status crpx_char_vector_new(size_t size, crpx_char_vector_t *out);
void crpx_char_vector_free(crpx_char_vector_t vector);
crpx_status crpx_char_vector_clone(const struct crpx_char_vector_struct *source, crpx_char_vector_t *out);
size_t crpx_char_vector_size(const struct crpx_char_vector_struct *vector);
const char *crpx_char_vector_get(const struct crpx_char_vector_struct *vector, size_t position);
size_t crpx_char_vector_length(const struct crpx_char_vector_struct *vector, size_t position);
crpx_status crpx_char_vector_resize(crpx_char_vector_t vector, size_t size);
crpx_status crpx_char_vector_set(crpx_char_vector_t vector, size_t position, const char *string);
crpx_status crpx_char_vector_append(crpx_char_vector_t vector, size_t position, const char *suffix);
crpx_status crpx_char_vector_push(crpx_char_vector_t vector, const char *string);
/* Reduce/reorder in one operation; repeated indices copy their strings. */
crpx_status crpx_char_vector_select(crpx_char_vector_t vector, const size_t *indices, size_t count);
/* Sort lexically or by decreasing length; old_indices maps new positions to old.
 * old_indices may be NULL. Ties retain original order. */
crpx_status crpx_char_vector_sort(crpx_char_vector_t vector, bool lexical, size_t *old_indices);
crpx_status crpx_char_vector_remove_empty(crpx_char_vector_t vector, size_t *removed);
crpx_status crpx_char_vector_deduplicate(crpx_char_vector_t vector, size_t *removed);
bool crpx_char_vector_equal(const struct crpx_char_vector_struct *a, const struct crpx_char_vector_struct *b);
crpx_status crpx_char_vector_find(const struct crpx_char_vector_struct *vector, const char *string, size_t *position);

#ifdef __cplusplus
}
#endif
#endif
