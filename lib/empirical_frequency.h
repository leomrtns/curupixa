/* Inspired by biomcmc-lib empirical_frequency.h.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_EMPIRICAL_FREQUENCY_H
#define CURUPIXA_EMPIRICAL_FREQUENCY_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_empfreq_struct *crpx_empfreq_t;
typedef struct {
  int64_t value;
  uint64_t frequency;
} crpx_empfreq_entry;

/* Objects own their storage and have no context or reference-count dependency.
 * Zero weights are omitted. Entries start ordered by decreasing frequency, with
 * ties resolved by increasing value. Empty histograms are valid. */
crpx_status crpx_empfreq_new(const int64_t *values, const uint64_t *weights, size_t size, crpx_empfreq_t *out);
crpx_status crpx_empfreq_merge(const struct crpx_empfreq_struct *a, const struct crpx_empfreq_struct *b,
                             crpx_empfreq_t *out);
void crpx_empfreq_free(crpx_empfreq_t frequencies);
size_t crpx_empfreq_size(const struct crpx_empfreq_struct *frequencies);
uint64_t crpx_empfreq_total(const struct crpx_empfreq_struct *frequencies);
/* Borrowed data remains valid until sort/free. Concurrent reads are safe. */
const crpx_empfreq_entry *crpx_empfreq_entries(const struct crpx_empfreq_struct *frequencies);
crpx_status crpx_empfreq_sort(crpx_empfreq_t frequencies, bool decreasing);
crpx_status crpx_empfreq_mode(const struct crpx_empfreq_struct *frequencies, int64_t *value);
crpx_status crpx_empfreq_range(const struct crpx_empfreq_struct *frequencies, int64_t *minimum, int64_t *maximum);
crpx_status crpx_mode_int64(const int64_t *values, const uint64_t *weights, size_t size, int64_t *out);

/* Return the original indices in value order, breaking ties by increasing index.
 * The output must have room for size indices and must not overlap the input. */
crpx_status crpx_rank_int64(const int64_t *values, size_t size, bool decreasing, size_t *indices);
crpx_status crpx_rank_double(const double *values, size_t size, bool decreasing, size_t *indices);

#ifdef __cplusplus
}
#endif
#endif
