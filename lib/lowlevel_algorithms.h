/* Reimplementation of biomcmc-lib lowlevel algorithms (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_LOWLEVEL_ALGORITHMS_H
#define CURUPIXA_LOWLEVEL_ALGORITHMS_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Rectangular minimum-cost assignment; rows <= columns. Finite, non-negative
 * row-major costs. Outputs one column per row and the total. Input is unchanged;
 * temporary storage O(rows+columns), time O(rows^2*columns). */
crpx_status crpx_assignment_solve(const double *costs, size_t rows, size_t columns,
                                 size_t *row_to_column, double *total);
/* Global edit distance, optionally allowing free trailing sequence ends.
 * Bytes are compared literally. Storage O(min(n,m)); arithmetic is checked. */
crpx_status crpx_levenshtein_distance(const void *a, size_t n, const void *b, size_t m,
                                     uint64_t substitution_cost, uint64_t indel_cost,
                                     bool free_ends, uint64_t *distance);
const char *crpx_string_last_substring(const char *text, const char *needle);
size_t crpx_string_common_prefix(const char *a, const char *b);

#ifdef __cplusplus
}
#endif
#endif
