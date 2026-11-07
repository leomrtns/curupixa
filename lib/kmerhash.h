/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_KMERHASH_H
#define CURUPIXA_KMERHASH_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_kmer_iterator crpx_kmer_iterator;
typedef struct {
  size_t position;
  uint64_t forward, reverse_complement, canonical, hash;
  bool reverse_is_canonical;
} crpx_kmer;
/* Exact two-bit encoding for k=1..32, with a seeded 64-bit final hash. A/C/G/T/U
 * are case-insensitive; every other byte breaks the window. Each iterator belongs to one worker. */
crpx_status crpx_kmer_iterator_new(size_t k, uint64_t seed, crpx_kmer_iterator **result);
void crpx_kmer_iterator_free(crpx_kmer_iterator *iterator);
/* Sequence bytes are borrowed and must remain valid and unchanged until reset/free. */
crpx_status crpx_kmer_iterator_reset(crpx_kmer_iterator *iterator, const char *sequence, size_t length);
/* NOT_FOUND signals exhaustion; output is unchanged on failure. */
crpx_status crpx_kmer_iterator_next(crpx_kmer_iterator *iterator, crpx_kmer *result);
#ifdef __cplusplus
}
#endif
#endif
