/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "rolling_hash.h"
#include <stdbool.h>
#include <stdlib.h>

struct crpx_rolling_hash {
  const unsigned char *data;
  size_t window, length, next;
  uint64_t hash, power, table[256];
  bool initialised;
};

#define CRPX_ROLLING_BASE UINT64_C(0x100000001b3)

static uint64_t power(uint64_t base, size_t exponent)
{
  uint64_t result = 1;
  while (exponent) {
    if (exponent & 1) result *= base;
    base *= base;
    exponent >>= 1;
  }
  return result;
}

static uint64_t next_seed(uint64_t *state)
{
  uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

crpx_status crpx_rolling_hash_new(size_t window, uint64_t seed, crpx_rolling_hash **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!window) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_rolling_hash *iterator = calloc(1, sizeof(*iterator));
  if (!iterator) return CRPX_ERROR_ALLOCATION;
  iterator->window = window;
  iterator->power = power(CRPX_ROLLING_BASE, window - 1);
  for (size_t i = 0; i < 256; i++) iterator->table[i] = next_seed(&seed);
  *result = iterator;
  return CRPX_OK;
}

void crpx_rolling_hash_free(crpx_rolling_hash *iterator) { free(iterator); }

crpx_status crpx_rolling_hash_reset(crpx_rolling_hash *iterator, const void *data, size_t length)
{
  if (!iterator || (!data && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  iterator->data = data;
  iterator->length = length;
  iterator->next = 0;
  iterator->hash = 0;
  iterator->initialised = false;
  return CRPX_OK;
}

crpx_status crpx_rolling_hash_next(crpx_rolling_hash *iterator, uint64_t *hash, size_t *position)
{
  if (!iterator || !hash || !position) return CRPX_ERROR_INVALID_ARGUMENT;
  if (iterator->window > iterator->length || iterator->next > iterator->length - iterator->window)
    return CRPX_ERROR_NOT_FOUND;
  if (!iterator->initialised) {
    for (size_t i = 0; i < iterator->window; i++)
      iterator->hash = iterator->hash * CRPX_ROLLING_BASE + iterator->table[iterator->data[i]];
    iterator->initialised = true;
  } else {
    unsigned char old = iterator->data[iterator->next - 1];
    unsigned char added = iterator->data[iterator->next + iterator->window - 1];
    iterator->hash = (iterator->hash - iterator->table[old] * iterator->power) * CRPX_ROLLING_BASE + iterator->table[added];
  }
  *hash = iterator->hash;
  *position = iterator->next++;
  return CRPX_OK;
}
