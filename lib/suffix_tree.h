/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_SUFFIX_TREE_H
#define CURUPIXA_SUFFIX_TREE_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif
/* A compact suffix-array/LCP index replaces the pointer-heavy suffix tree. Construction copies
 * arbitrary bytes, uses O(n) memory and O(n log^2 n) time. Queries may run concurrently. */
typedef struct crpx_suffix_index crpx_suffix_index;
crpx_status crpx_suffix_index_new(const void *data, size_t length, crpx_suffix_index **result);
void crpx_suffix_index_free(crpx_suffix_index *index);
size_t crpx_suffix_index_length(const crpx_suffix_index *index);
/* Find a nonempty byte pattern, returning a contiguous suffix-array range. NOT_FOUND leaves
 * outputs unchanged. Match positions follow lexical suffix order, rather than coordinate order. */
crpx_status crpx_suffix_index_find(const crpx_suffix_index *index, const void *pattern, size_t length,
                                 size_t *first, size_t *count);
crpx_status crpx_suffix_index_position(const crpx_suffix_index *index, size_t suffix_rank, size_t *position);
/* Longest common prefix of the suffixes at two sequence coordinates (coordinate length is empty). */
crpx_status crpx_suffix_index_lcp(const crpx_suffix_index *index, size_t first_position,
                                size_t second_position, size_t *length);
/* Longest repeated substring, allowing overlapping occurrences; zero length for no repeats. */
crpx_status crpx_suffix_index_longest_repeat(const crpx_suffix_index *index, size_t *position, size_t *length);
#ifdef __cplusplus
}
#endif
#endif
