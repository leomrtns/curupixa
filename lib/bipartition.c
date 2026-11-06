/* Inspired by biomcmc-lib bipartition.c (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "bipartition.h"
#include <stdlib.h>
#include <string.h>

struct crpx_bipartition_struct {
  size_t size, nwords;
  uint64_t *words;
};

static uint64_t last_mask(size_t size)
{
  return size % 64 ? (UINT64_C(1) << (size % 64)) - 1 : UINT64_MAX;
}

crpx_status crpx_bipartition_new(size_t size, crpx_bipartition_t *out)
{
  crpx_bipartition_t b;
  size_t count = size / 64 + (size % 64 != 0);
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (count > SIZE_MAX / sizeof(uint64_t)) return CRPX_ERROR_OVERFLOW;
  b = calloc(1, sizeof(*b));
  if (!b) return CRPX_ERROR_ALLOCATION;
  if (count && !(b->words = calloc(count, sizeof(*b->words)))) { free(b); return CRPX_ERROR_ALLOCATION; }
  b->size = size;
  b->nwords = count;
  *out = b;
  return CRPX_OK;
}

crpx_status crpx_bipartition_clone(const struct crpx_bipartition_struct *source, crpx_bipartition_t *out)
{
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!source) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_bipartition_new(source->size, out);
  if (status == CRPX_OK && source->nwords) memcpy((*out)->words, source->words, source->nwords * sizeof(uint64_t));
  return status;
}

void crpx_bipartition_free(crpx_bipartition_t bits)
{
  if (bits) { free(bits->words); free(bits); }
}

size_t crpx_bipartition_size(const struct crpx_bipartition_struct *bits) { return bits ? bits->size : 0; }

size_t crpx_bipartition_count(const struct crpx_bipartition_struct *bits)
{
  size_t result = 0;
  if (bits) for (size_t i = 0; i < bits->nwords; ++i) {
    uint64_t value = bits->words[i];
    /* Portable SWAR popcount: compilers may replace this with a CPU instruction. */
    value -= (value >> 1) & UINT64_C(0x5555555555555555);
    value = (value & UINT64_C(0x3333333333333333)) + ((value >> 2) & UINT64_C(0x3333333333333333));
    value = (value + (value >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
    result += (value * UINT64_C(0x0101010101010101)) >> 56;
  }
  return result;
}

crpx_status crpx_bipartition_resize(crpx_bipartition_t bits, size_t size)
{
  size_t count = size / 64 + (size % 64 != 0);
  uint64_t *words;
  if (!bits) return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / sizeof(uint64_t)) return CRPX_ERROR_OVERFLOW;
  if (count != bits->nwords) {
    words = count ? calloc(count, sizeof(*words)) : NULL;
    if (count && !words) return CRPX_ERROR_ALLOCATION;
    size_t common = count < bits->nwords ? count : bits->nwords;
    if (common) memcpy(words, bits->words, common * sizeof(*words));
    free(bits->words);
    bits->words = words;
    bits->nwords = count;
  }
  bits->size = size;
  if (count) bits->words[count - 1] &= last_mask(size);
  return CRPX_OK;
}

crpx_status crpx_bipartition_set(crpx_bipartition_t bits, size_t position, bool value)
{
  if (!bits) return CRPX_ERROR_INVALID_ARGUMENT;
  if (position >= bits->size) return CRPX_ERROR_OUT_OF_RANGE;
  uint64_t mask = UINT64_C(1) << (position % 64);
  if (value) bits->words[position / 64] |= mask;
  else bits->words[position / 64] &= ~mask;
  return CRPX_OK;
}

crpx_status crpx_bipartition_get(const struct crpx_bipartition_struct *bits, size_t position, bool *value)
{
  if (!bits || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  if (position >= bits->size) return CRPX_ERROR_OUT_OF_RANGE;
  *value = (bits->words[position / 64] & (UINT64_C(1) << (position % 64))) != 0;
  return CRPX_OK;
}

void crpx_bipartition_clear(crpx_bipartition_t bits)
{
  if (bits && bits->nwords) memset(bits->words, 0, bits->nwords * sizeof(uint64_t));
}

crpx_status crpx_bipartition_singleton(crpx_bipartition_t bits, size_t position)
{
  if (!bits) return CRPX_ERROR_INVALID_ARGUMENT;
  if (position >= bits->size) return CRPX_ERROR_OUT_OF_RANGE;
  crpx_bipartition_clear(bits);
  return crpx_bipartition_set(bits, position, true);
}

crpx_status crpx_bipartition_copy(crpx_bipartition_t to, const struct crpx_bipartition_struct *from)
{
  if (!to || !from || to->size != from->size) return CRPX_ERROR_INVALID_ARGUMENT;
  if (to != from && to->nwords) memcpy(to->words, from->words, to->nwords * sizeof(uint64_t));
  return CRPX_OK;
}

crpx_status crpx_bipartition_combine(crpx_bipartition_t out, const struct crpx_bipartition_struct *a,
                                   const struct crpx_bipartition_struct *b, crpx_bit_operation op)
{
  if (!out || !a || !b || out->size != a->size || a->size != b->size ||
      op < CRPX_BITS_OR || op > CRPX_BITS_XORNOT) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < out->nwords; ++i) {
    uint64_t x = a->words[i], y = b->words[i];
    switch (op) {
      case CRPX_BITS_OR: out->words[i] = x | y; break;
      case CRPX_BITS_AND: out->words[i] = x & y; break;
      case CRPX_BITS_ANDNOT: out->words[i] = x & ~y; break;
      case CRPX_BITS_XOR: out->words[i] = x ^ y; break;
      case CRPX_BITS_XORNOT: out->words[i] = x ^ ~y; break;
    }
  }
  if (out->nwords) out->words[out->nwords - 1] &= last_mask(out->size);
  return CRPX_OK;
}

crpx_status crpx_bipartition_complement(crpx_bipartition_t out, const struct crpx_bipartition_struct *source)
{
  if (!out || !source || out->size != source->size) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < out->nwords; ++i) out->words[i] = ~source->words[i];
  if (out->nwords) out->words[out->nwords - 1] &= last_mask(out->size);
  return CRPX_OK;
}

bool crpx_bipartition_equal(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b)
{
  return a && b && a->size == b->size && (!a->nwords || !memcmp(a->words, b->words, a->nwords * sizeof(uint64_t)));
}

bool crpx_bipartition_same_split(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b)
{
  if (!a || !b || a->size != b->size) return false;
  if (crpx_bipartition_equal(a, b)) return true;
  for (size_t i = 0; i < a->nwords; ++i) {
    uint64_t mask = i + 1 == a->nwords ? last_mask(a->size) : UINT64_MAX;
    if ((a->words[i] ^ b->words[i]) != mask) return false;
  }
  return true;
}

bool crpx_bipartition_contains(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b)
{
  if (!a || !b || a->size != b->size) return false;
  for (size_t i = 0; i < a->nwords; ++i) if ((a->words[i] & b->words[i]) != b->words[i]) return false;
  return true;
}

bool crpx_bipartition_compatible(const struct crpx_bipartition_struct *a, const struct crpx_bipartition_struct *b)
{
  bool nonempty[4] = {false, false, false, false};
  if (!a || !b || a->size != b->size) return false;
  for (size_t i = 0; i < a->nwords; ++i) {
    uint64_t x = a->words[i], y = b->words[i];
    uint64_t mask = i + 1 == a->nwords ? last_mask(a->size) : UINT64_MAX;
    nonempty[0] |= (x & y) != 0;
    nonempty[1] |= (x & ~y & mask) != 0;
    nonempty[2] |= (~x & y & mask) != 0;
    nonempty[3] |= (~x & ~y & mask) != 0;
  }
  return !nonempty[0] || !nonempty[1] || !nonempty[2] || !nonempty[3];
}

void crpx_bipartition_canonicalise(crpx_bipartition_t bits)
{
  if (!bits || !bits->size) return;
  size_t count = crpx_bipartition_count(bits);
  if (count > bits->size - count || (count == bits->size - count && (bits->words[0] & 1)))
    (void) crpx_bipartition_complement(bits, bits);
}

crpx_status crpx_bipartition_indices(const struct crpx_bipartition_struct *bits, size_t *indices,
                                   size_t capacity, size_t *required)
{
  size_t count;
  if (!bits || !required || (!indices && capacity)) return CRPX_ERROR_INVALID_ARGUMENT;
  *required = count = crpx_bipartition_count(bits);
  if (capacity < count) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t i = 0, j = 0; i < bits->size; ++i)
    if (bits->words[i / 64] & (UINT64_C(1) << (i % 64))) indices[j++] = i;
  return CRPX_OK;
}

const uint64_t *crpx_bipartition_words(const struct crpx_bipartition_struct *bits, size_t *count)
{
  if (count) *count = bits ? bits->nwords : 0;
  return bits ? bits->words : NULL;
}
