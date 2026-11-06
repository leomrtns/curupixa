/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Classic HyperLogLog estimator from Flajolet, Fusy, Gandouet and Meunier (2007):
 * https://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf
 * Fresh implementation; no mutable lookup tables. */
#include "hll.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_hll {
  uint64_t seed;
  unsigned precision;
  size_t count;
  uint8_t registers[];
};

static uint64_t mix(uint64_t value)
{
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

crpx_status crpx_hll_new(unsigned precision, uint64_t seed, crpx_hll **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (precision < 4 || precision > 18) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t count = (size_t)1 << precision;
  crpx_hll *sketch = calloc(1, sizeof(*sketch) + count);
  if (!sketch) return CRPX_ERROR_ALLOCATION;
  sketch->precision = precision;
  sketch->seed = seed;
  sketch->count = count;
  *result = sketch;
  return CRPX_OK;
}

void crpx_hll_free(crpx_hll *sketch) { free(sketch); }
void crpx_hll_reset(crpx_hll *sketch) { if (sketch) memset(sketch->registers, 0, sketch->count); }
size_t crpx_hll_register_count(const crpx_hll *sketch) { return sketch ? sketch->count : 0; }

crpx_status crpx_hll_clone(const crpx_hll *source, crpx_hll **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!source) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_hll *copy;
  crpx_status status = crpx_hll_new(source->precision, source->seed, &copy);
  if (status != CRPX_OK) return status;
  memcpy(copy->registers, source->registers, source->count);
  *result = copy;
  return CRPX_OK;
}

crpx_status crpx_hll_add_hash(crpx_hll *sketch, uint64_t hash)
{
  if (!sketch) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t bucket = hash & (sketch->count - 1);
  uint64_t tail = hash >> sketch->precision;
  uint8_t rank = 1;
  /* The bucket occupies the low bits, so count trailing zeroes of the remaining bits. */
  while (rank <= 64 - sketch->precision && !(tail & 1)) { rank++; tail >>= 1; }
  if (rank > sketch->registers[bucket]) sketch->registers[bucket] = rank;
  return CRPX_OK;
}

crpx_status crpx_hll_add(crpx_hll *sketch, const void *data, size_t length)
{
  if (!sketch || (!data && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  const unsigned char *bytes = data;
  uint64_t hash = UINT64_C(0xcbf29ce484222325) ^ sketch->seed;
  for (size_t i = 0; i < length; i++) hash = (hash ^ bytes[i]) * UINT64_C(0x100000001b3);
  return crpx_hll_add_hash(sketch, mix(hash));
}

crpx_status crpx_hll_merge(crpx_hll *destination, const crpx_hll *source)
{
  if (!destination || !source) return CRPX_ERROR_INVALID_ARGUMENT;
  if (destination->precision != source->precision || destination->seed != source->seed)
    return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < destination->count; i++)
    if (source->registers[i] > destination->registers[i]) destination->registers[i] = source->registers[i];
  return CRPX_OK;
}

crpx_status crpx_hll_estimate(const crpx_hll *sketch, double *estimate, double *relative_standard_error)
{
  if (!sketch || !estimate) return CRPX_ERROR_INVALID_ARGUMENT;
  double count = (double)sketch->count;
  double alpha = 0.7213 / (1 + 1.079 / count);
  if (sketch->count == 16) alpha = 0.673;
  else if (sketch->count == 32) alpha = 0.697;
  else if (sketch->count == 64) alpha = 0.709;
  double sum = 0;
  size_t empty = 0;
  for (size_t i = 0; i < sketch->count; i++) {
    sum += ldexp(1.0, -(int)sketch->registers[i]);
    empty += sketch->registers[i] == 0;
  }
  double value = alpha * count * count / sum;
  if (value <= 2.5 * count && empty) value = count * log(count / (double)empty);
  else if (value > 0x1.0p64 / 30) {
    if (value >= 0x1.0p64) return CRPX_ERROR_OVERFLOW;
    value = -0x1.0p64 * log1p(-value / 0x1.0p64);
  }
  *estimate = value;
  if (relative_standard_error) *relative_standard_error = 1.04 / sqrt(count);
  return CRPX_OK;
}

crpx_status crpx_hll_serialize(const crpx_hll *sketch, void *buffer, size_t capacity, size_t *written)
{
  if (!sketch || !written) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t required = 24 + sketch->count;
  *written = required;
  if (!buffer || capacity < required) return CRPX_ERROR_OUT_OF_RANGE;
  unsigned char *bytes = buffer;
  memset(bytes, 0, 24);
  memcpy(bytes, "CRPXHLL1", 8);
  bytes[8] = (unsigned char)sketch->precision;
  for (size_t i = 0; i < 8; i++) bytes[16 + i] = (unsigned char)(sketch->seed >> (8 * i));
  memcpy(bytes + 24, sketch->registers, sketch->count);
  return CRPX_OK;
}

crpx_status crpx_hll_deserialize(const void *buffer, size_t length, crpx_hll **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!buffer) return CRPX_ERROR_INVALID_ARGUMENT;
  const unsigned char *bytes = buffer;
  if (length < 24 || memcmp(bytes, "CRPXHLL1", 8)) return CRPX_ERROR_PARSE;
  unsigned precision = bytes[8];
  if (precision < 4 || precision > 18) return CRPX_ERROR_PARSE;
  size_t count = (size_t)1 << precision;
  if (length != 24 + count) return CRPX_ERROR_PARSE;
  for (size_t i = 9; i < 16; i++) if (bytes[i]) return CRPX_ERROR_PARSE;
  for (size_t i = 0; i < count; i++) if (bytes[24 + i] > 65 - precision) return CRPX_ERROR_PARSE;
  uint64_t seed = 0;
  for (size_t i = 0; i < 8; i++) seed |= (uint64_t)bytes[16 + i] << (8 * i);
  crpx_hll *sketch;
  crpx_status status = crpx_hll_new(precision, seed, &sketch);
  if (status != CRPX_OK) return status;
  memcpy(sketch->registers, bytes + 24, count);
  *result = sketch;
  return CRPX_OK;
}
