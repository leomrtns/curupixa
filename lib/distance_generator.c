/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "distance_generator.h"
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct crpx_distance_generator_struct {
  size_t samples, components, pairs;
  crpx_distance_callback callback;
  void *data;
  double **cache;
  atomic_flag lock;
};

static void lock_generator(crpx_distance_generator_t generator)
{ while (atomic_flag_test_and_set_explicit(&generator->lock, memory_order_acquire)) { } }
static void unlock_generator(crpx_distance_generator_t generator)
{ atomic_flag_clear_explicit(&generator->lock, memory_order_release); }
static void reset_generator(crpx_distance_generator_t generator)
{
  if (generator->cache)
    for (size_t i = 0; i < generator->pairs; ++i) { free(generator->cache[i]); generator->cache[i] = NULL; }
}

crpx_status crpx_distance_generator_new(size_t samples, size_t components, crpx_distance_callback callback,
                                       void *data, bool cache, crpx_distance_generator_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!samples || !components || !callback) return CRPX_ERROR_INVALID_ARGUMENT;
  if (components > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  size_t first = samples, second = samples - 1;
  if (first % 2) second /= 2; else first /= 2;
  if (first && second > SIZE_MAX / first) return CRPX_ERROR_OVERFLOW;
  size_t pairs = first * second;
  if (cache && pairs > SIZE_MAX / sizeof(double *)) return CRPX_ERROR_OVERFLOW;
  crpx_distance_generator_t generator = calloc(1, sizeof(*generator));
  if (!generator) return CRPX_ERROR_ALLOCATION;
  generator->samples = samples; generator->components = components; generator->pairs = pairs;
  generator->callback = callback; generator->data = data;
  atomic_flag_clear(&generator->lock);
  if (cache && pairs) {
    generator->cache = calloc(pairs, sizeof(double *));
    if (!generator->cache) { free(generator); return CRPX_ERROR_ALLOCATION; }
  }
  *out = generator;
  return CRPX_OK;
}

void crpx_distance_generator_free(crpx_distance_generator_t generator)
{
  if (generator) { reset_generator(generator); free(generator->cache); free(generator); }
}
size_t crpx_distance_generator_samples(const struct crpx_distance_generator_struct *generator)
{ return generator ? generator->samples : 0; }
size_t crpx_distance_generator_components(const struct crpx_distance_generator_struct *generator)
{ return generator ? generator->components : 0; }

static crpx_status calculate_locked(crpx_distance_generator_t generator, size_t first, size_t second,
                                     double **values, bool *temporary)
{
  if (first > second) { size_t swap = first; first = second; second = swap; }
  size_t a = second, b = second - 1;
  if (a % 2) b /= 2; else a /= 2;
  size_t index = a * b + first;
  *temporary = false;
  if (generator->cache && generator->cache[index]) { *values = generator->cache[index]; return CRPX_OK; }
  double *result = malloc(generator->components * sizeof(double));
  if (!result) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < generator->components; ++i) result[i] = NAN;
  crpx_status status = generator->callback(generator->data, first, second, result, generator->components);
  if (status == CRPX_OK)
    for (size_t i = 0; i < generator->components; ++i)
      if (!isfinite(result[i])) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
  if (status != CRPX_OK) { free(result); return status; }
  if (generator->cache) generator->cache[index] = result;
  else *temporary = true;
  *values = result;
  return CRPX_OK;
}

crpx_status crpx_distance_generator_get(crpx_distance_generator_t generator, size_t first, size_t second,
                                       size_t component, double *out)
{
  if (!generator || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (first >= generator->samples || second >= generator->samples || component >= generator->components)
    return CRPX_ERROR_OUT_OF_RANGE;
  if (first == second) { *out = 0; return CRPX_OK; }
  lock_generator(generator);
  double *values; bool temporary;
  crpx_status status = calculate_locked(generator, first, second, &values, &temporary);
  if (status == CRPX_OK) { *out = values[component]; if (temporary) free(values); }
  unlock_generator(generator);
  return status;
}

crpx_status crpx_distance_generator_get_all(crpx_distance_generator_t generator, size_t first, size_t second,
                                           double *out, size_t capacity)
{
  if (!generator || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (first >= generator->samples || second >= generator->samples || capacity < generator->components)
    return CRPX_ERROR_OUT_OF_RANGE;
  if (first == second) { for (size_t i = 0; i < generator->components; ++i) out[i] = 0; return CRPX_OK; }
  lock_generator(generator);
  double *values; bool temporary;
  crpx_status status = calculate_locked(generator, first, second, &values, &temporary);
  if (status == CRPX_OK) { memcpy(out, values, generator->components * sizeof(double)); if (temporary) free(values); }
  unlock_generator(generator);
  return status;
}

crpx_status crpx_distance_generator_reset(crpx_distance_generator_t generator)
{
  if (!generator) return CRPX_ERROR_INVALID_ARGUMENT;
  lock_generator(generator); reset_generator(generator); unlock_generator(generator);
  return CRPX_OK;
}

crpx_status crpx_distance_generator_set_callback(crpx_distance_generator_t generator,
                                                crpx_distance_callback callback, void *data)
{
  if (!generator || !callback) return CRPX_ERROR_INVALID_ARGUMENT;
  lock_generator(generator);
  reset_generator(generator);
  generator->callback = callback; generator->data = data;
  unlock_generator(generator);
  return CRPX_OK;
}
