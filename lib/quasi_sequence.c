/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "quasi_sequence.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
struct crpx_quasi_sequence {
  crpx_quasi_kind kind;
  size_t dimensions;
  uint64_t index, *directions;
  unsigned *bases;
};
static bool is_prime (unsigned number);

crpx_status
crpx_quasi_sequence_new (crpx_quasi_kind kind, size_t dimensions, crpx_quasi_sequence **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!dimensions || kind < CRPX_QUASI_HALTON || kind > CRPX_QUASI_SOBOL) return CRPX_ERROR_INVALID_ARGUMENT;
  if (dimensions > (kind == CRPX_QUASI_SOBOL ? 8u : 2048u)) return CRPX_ERROR_UNSUPPORTED;
  crpx_quasi_sequence *sequence = calloc (1, sizeof (*sequence));
  if (!sequence) return CRPX_ERROR_ALLOCATION;
  sequence->kind = kind;
  sequence->dimensions = dimensions;
  sequence->index = 1;
  if (kind == CRPX_QUASI_SOBOL) {
    static const unsigned degree[] = {0, 1, 2, 3, 3, 4, 4, 5};
    static const unsigned polynomial[] = {0, 0, 1, 1, 2, 1, 4, 2};
    static const unsigned initial[8][5] = {{0},       {1},          {1, 3},        {1, 3, 1},
                                           {1, 1, 1}, {1, 1, 3, 3}, {1, 3, 5, 13}, {1, 1, 5, 5, 17}};
    sequence->directions = calloc (dimensions * 64, sizeof (uint64_t));
    if (!sequence->directions) {
      crpx_quasi_sequence_free (sequence);
      return CRPX_ERROR_ALLOCATION;
    }
    for (size_t dimension = 0; dimension < dimensions; dimension++) {
      uint64_t *directions = sequence->directions + dimension * 64;
      if (!dimension)
        for (unsigned bit = 0; bit < 64; bit++)
          directions[bit] = UINT64_C (1) << (63 - bit);
      else {
        unsigned order = degree[dimension];
        for (unsigned bit = 0; bit < order; bit++)
          directions[bit] = (uint64_t) initial[dimension][bit] << (63 - bit);
        for (unsigned bit = order; bit < 64; bit++) {
          directions[bit] = directions[bit - order] ^ (directions[bit - order] >> order);
          for (unsigned term = 1; term < order; term++)
            if ((polynomial[dimension] >> (order - 1 - term)) & 1) directions[bit] ^= directions[bit - term];
        }
      }
    }
  } else {
    sequence->bases = malloc (dimensions * sizeof (*sequence->bases));
    if (!sequence->bases) {
      crpx_quasi_sequence_free (sequence);
      return CRPX_ERROR_ALLOCATION;
    }
    unsigned prime = 2;
    for (size_t i = 0; i < dimensions; i++) {
      while (!is_prime (prime))
        prime++;
      sequence->bases[i] = prime++;
    }
  }
  *out = sequence;
  return CRPX_OK;
}
void
crpx_quasi_sequence_free (crpx_quasi_sequence *sequence)
{
  if (sequence) {
    free (sequence->directions);
    free (sequence->bases);
    free (sequence);
  }
}
crpx_status
crpx_quasi_sequence_seek (crpx_quasi_sequence *sequence, uint64_t next_index)
{
  if (!sequence || next_index >= UINT64_C (9007199254740992)) return CRPX_ERROR_INVALID_ARGUMENT;
  sequence->index = next_index;
  return CRPX_OK;
}
crpx_status
crpx_quasi_sequence_next (crpx_quasi_sequence *sequence, double *point, size_t capacity)
{ /* Seek evaluates a point directly, so independent workers can request disjoint index ranges. */
  if (!sequence || !point || capacity < sequence->dimensions) return CRPX_ERROR_INVALID_ARGUMENT;
  if (sequence->index >= UINT64_C (9007199254740992)) return CRPX_ERROR_OVERFLOW;
  for (size_t dimension = 0; dimension < sequence->dimensions; dimension++) {
    uint64_t index = sequence->index;
    if (sequence->kind == CRPX_QUASI_SOBOL) {
      uint64_t gray = index ^ (index >> 1), value = 0;
      for (unsigned bit = 0; gray; bit++, gray >>= 1)
        if (gray & 1) value ^= sequence->directions[dimension * 64 + bit];
      point[dimension] = (double) (value >> 11) * 0x1p-53;
    } else {
      unsigned base = sequence->bases[dimension];
      long double value = 0, place = 1.0L / base;
      while (index) {
        unsigned digit = (unsigned) (index % base);
        if (sequence->kind == CRPX_QUASI_REVERSE_HALTON && digit) digit = base - digit;
        value += digit * place;
        place /= base;
        index /= base;
      }
      point[dimension] = fmin ((double) value, nextafter (1.0, 0.0));
    }
  }
  sequence->index++;
  return CRPX_OK;
}
static bool
is_prime (unsigned number)
{
  for (unsigned divisor = 2; divisor <= number / divisor; divisor++)
    if (!(number % divisor)) return false;
  return number >= 2;
}
