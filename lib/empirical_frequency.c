/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "empirical_frequency.h"
#include <math.h>
#include <stdlib.h>

struct crpx_empfreq_struct {
  size_t size;
  uint64_t total;
  crpx_empfreq_entry *entries;
};

static int compare_value(const void *left, const void *right)
{
  const crpx_empfreq_entry *a = left, *b = right;
  return (a->value > b->value) - (a->value < b->value);
}

static int compare_frequency_ascending(const void *left, const void *right)
{
  const crpx_empfreq_entry *a = left, *b = right;
  if (a->frequency != b->frequency) return (a->frequency > b->frequency) ? 1 : -1;
  return compare_value(left, right);
}

static int compare_frequency_descending(const void *left, const void *right)
{
  const crpx_empfreq_entry *a = left, *b = right;
  if (a->frequency != b->frequency) return (a->frequency < b->frequency) ? 1 : -1;
  return compare_value(left, right);
}

static crpx_status new_entries(size_t capacity, crpx_empfreq_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (capacity > SIZE_MAX / sizeof(crpx_empfreq_entry)) return CRPX_ERROR_OVERFLOW;
  crpx_empfreq_t frequencies = calloc(1, sizeof(*frequencies));
  if (!frequencies) return CRPX_ERROR_ALLOCATION;
  if (capacity) {
    frequencies->entries = malloc(capacity * sizeof(*frequencies->entries));
    if (!frequencies->entries) { free(frequencies); return CRPX_ERROR_ALLOCATION; }
  }
  *out = frequencies;
  return CRPX_OK;
}

static crpx_status aggregate(crpx_empfreq_t frequencies)
{
  if (frequencies->size > 1)
    qsort(frequencies->entries, frequencies->size, sizeof(*frequencies->entries), compare_value);
  size_t count = 0;
  uint64_t total = 0;
  for (size_t i = 0; i < frequencies->size; ++i) {
    crpx_empfreq_entry entry = frequencies->entries[i];
    if (entry.frequency > UINT64_MAX - total) return CRPX_ERROR_OVERFLOW;
    total += entry.frequency;
    if (count && frequencies->entries[count - 1].value == entry.value)
      frequencies->entries[count - 1].frequency += entry.frequency;
    else frequencies->entries[count++] = entry;
  }
  frequencies->size = count;
  frequencies->total = total;
  return crpx_empfreq_sort(frequencies, true);
}

crpx_status crpx_empfreq_new(const int64_t *values, const uint64_t *weights, size_t size, crpx_empfreq_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (size && !values) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_empfreq_t frequencies;
  crpx_status status = new_entries(size, &frequencies);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < size; ++i) {
    uint64_t weight = weights ? weights[i] : 1;
    if (weight) frequencies->entries[frequencies->size++] = (crpx_empfreq_entry){values[i], weight};
  }
  status = aggregate(frequencies);
  if (status != CRPX_OK) { crpx_empfreq_free(frequencies); return status; }
  *out = frequencies;
  return CRPX_OK;
}

crpx_status crpx_empfreq_merge(const struct crpx_empfreq_struct *a, const struct crpx_empfreq_struct *b,
                             crpx_empfreq_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!a || !b) return CRPX_ERROR_INVALID_ARGUMENT;
  if (a->size > SIZE_MAX - b->size) return CRPX_ERROR_OVERFLOW;
  crpx_empfreq_t frequencies;
  crpx_status status = new_entries(a->size + b->size, &frequencies);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < a->size; ++i) frequencies->entries[frequencies->size++] = a->entries[i];
  for (size_t i = 0; i < b->size; ++i) frequencies->entries[frequencies->size++] = b->entries[i];
  status = aggregate(frequencies);
  if (status != CRPX_OK) { crpx_empfreq_free(frequencies); return status; }
  *out = frequencies;
  return CRPX_OK;
}

void crpx_empfreq_free(crpx_empfreq_t frequencies)
{
  if (frequencies) { free(frequencies->entries); free(frequencies); }
}

size_t crpx_empfreq_size(const struct crpx_empfreq_struct *frequencies)
{ return frequencies ? frequencies->size : 0; }

uint64_t crpx_empfreq_total(const struct crpx_empfreq_struct *frequencies)
{ return frequencies ? frequencies->total : 0; }

const crpx_empfreq_entry *crpx_empfreq_entries(const struct crpx_empfreq_struct *frequencies)
{ return frequencies ? frequencies->entries : NULL; }

crpx_status crpx_empfreq_sort(crpx_empfreq_t frequencies, bool decreasing)
{
  if (!frequencies) return CRPX_ERROR_INVALID_ARGUMENT;
  if (frequencies->size > 1)
    qsort(frequencies->entries, frequencies->size, sizeof(*frequencies->entries),
          decreasing ? compare_frequency_descending : compare_frequency_ascending);
  return CRPX_OK;
}

crpx_status crpx_empfreq_mode(const struct crpx_empfreq_struct *frequencies, int64_t *value)
{
  if (!frequencies || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!frequencies->size) return CRPX_ERROR_NOT_FOUND;
  const crpx_empfreq_entry *best = frequencies->entries;
  for (size_t i = 1; i < frequencies->size; ++i)
    if (compare_frequency_descending(frequencies->entries + i, best) < 0) best = frequencies->entries + i;
  *value = best->value;
  return CRPX_OK;
}

crpx_status crpx_empfreq_range(const struct crpx_empfreq_struct *frequencies, int64_t *minimum, int64_t *maximum)
{
  if (!frequencies || !minimum || !maximum) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!frequencies->size) return CRPX_ERROR_NOT_FOUND;
  int64_t low = frequencies->entries[0].value, high = low;
  for (size_t i = 1; i < frequencies->size; ++i) {
    if (frequencies->entries[i].value < low) low = frequencies->entries[i].value;
    if (frequencies->entries[i].value > high) high = frequencies->entries[i].value;
  }
  *minimum = low;
  *maximum = high;
  return CRPX_OK;
}

crpx_status crpx_mode_int64(const int64_t *values, const uint64_t *weights, size_t size, int64_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_empfreq_t frequencies;
  crpx_status status = crpx_empfreq_new(values, weights, size, &frequencies);
  if (status != CRPX_OK) return status;
  status = crpx_empfreq_mode(frequencies, out);
  crpx_empfreq_free(frequencies);
  return status;
}

typedef struct { int64_t value; size_t index; } ranked_int64;
typedef struct { double value; size_t index; } ranked_double;

#define DEFINE_RANK(TYPE, NAME) \
static int compare_##NAME##_ascending(const void *left, const void *right) \
{ \
  const ranked_##NAME *a = left, *b = right; \
  if (a->value != b->value) return (a->value > b->value) ? 1 : -1; \
  return (a->index > b->index) - (a->index < b->index); \
} \
static int compare_##NAME##_descending(const void *left, const void *right) \
{ \
  const ranked_##NAME *a = left, *b = right; \
  if (a->value != b->value) return (a->value < b->value) ? 1 : -1; \
  return (a->index > b->index) - (a->index < b->index); \
} \
crpx_status crpx_rank_##NAME(const TYPE *values, size_t size, bool decreasing, size_t *indices) \
{ \
  if (size && (!values || !indices)) return CRPX_ERROR_INVALID_ARGUMENT; \
  if (size > SIZE_MAX / sizeof(ranked_##NAME)) return CRPX_ERROR_OVERFLOW; \
  if (!size) return CRPX_OK; \
  ranked_##NAME *ranked = malloc(size * sizeof(*ranked)); \
  if (!ranked) return CRPX_ERROR_ALLOCATION; \
  for (size_t i = 0; i < size; ++i) { \
    if (isnan((double)values[i])) { free(ranked); return CRPX_ERROR_INVALID_ARGUMENT; } \
    ranked[i] = (ranked_##NAME){values[i], i}; \
  } \
  qsort(ranked, size, sizeof(*ranked), decreasing ? compare_##NAME##_descending : compare_##NAME##_ascending); \
  for (size_t i = 0; i < size; ++i) indices[i] = ranked[i].index; \
  free(ranked); \
  return CRPX_OK; \
}

DEFINE_RANK(int64_t, int64)
DEFINE_RANK(double, double)
