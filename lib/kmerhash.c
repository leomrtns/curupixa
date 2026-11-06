/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "kmerhash.h"
#include <stdlib.h>

struct crpx_kmer_iterator {
  const unsigned char *sequence;
  size_t length, position, valid, k;
  uint64_t seed, forward, reverse, mask;
};

static unsigned base_code(unsigned char base)
{
  switch (base) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': case 'U': case 'u': return 3;
    default: return 4;
  }
}

static uint64_t mix(uint64_t value)
{
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

crpx_status crpx_kmer_iterator_new(size_t k, uint64_t seed, crpx_kmer_iterator **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (k < 1 || k > 32) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_kmer_iterator *iterator = calloc(1, sizeof(*iterator));
  if (!iterator) return CRPX_ERROR_ALLOCATION;
  iterator->k = k;
  iterator->seed = seed;
  iterator->mask = k == 32 ? UINT64_MAX : (UINT64_C(1) << (2 * k)) - 1;
  *result = iterator;
  return CRPX_OK;
}

void crpx_kmer_iterator_free(crpx_kmer_iterator *iterator) { free(iterator); }

crpx_status crpx_kmer_iterator_reset(crpx_kmer_iterator *iterator, const char *sequence, size_t length)
{
  if (!iterator || (!sequence && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  iterator->sequence = (const unsigned char *)sequence;
  iterator->length = length;
  iterator->position = iterator->valid = 0;
  iterator->forward = iterator->reverse = 0;
  return CRPX_OK;
}

crpx_status crpx_kmer_iterator_next(crpx_kmer_iterator *iterator, crpx_kmer *result)
{
  if (!iterator || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  while (iterator->position < iterator->length) {
    unsigned code = base_code(iterator->sequence[iterator->position++]);
    if (code == 4) {
      iterator->valid = 0;
      iterator->forward = iterator->reverse = 0;
      continue;
    }
    iterator->forward = ((iterator->forward << 2) | code) & iterator->mask;
    iterator->reverse = (iterator->reverse >> 2) | (uint64_t)(3 - code) << (2 * (iterator->k - 1));
    if (iterator->valid < iterator->k) iterator->valid++;
    if (iterator->valid < iterator->k) continue;
    crpx_kmer value;
    value.position = iterator->position - iterator->k;
    value.forward = iterator->forward;
    value.reverse_complement = iterator->reverse;
    value.reverse_is_canonical = value.reverse_complement < value.forward;
    value.canonical = value.reverse_is_canonical ? value.reverse_complement : value.forward;
    value.hash = mix(mix(value.canonical ^ iterator->seed) ^ (iterator->k * UINT64_C(0x9e3779b97f4a7c15)));
    *result = value;
    return CRPX_OK;
  }
  return CRPX_ERROR_NOT_FOUND;
}
