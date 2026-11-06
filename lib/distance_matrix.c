/* Reimplementation of biomcmc-lib distance_matrix (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "distance_matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_distance_matrix_struct { size_t size, count; double *values; };
struct species_pair { double mean, minimum; uint64_t count; };
struct crpx_species_distances_struct { size_t size, count; struct species_pair *pairs; };

static crpx_status pair_count(size_t size, size_t *count)
{
  size_t a = size, b = size ? size - 1 : 0;
  if (a % 2) b /= 2; else a /= 2;
  if (b && a > SIZE_MAX / b) return CRPX_ERROR_OVERFLOW;
  *count = a * b;
  return CRPX_OK;
}
static size_t pair_index(size_t a, size_t b)
{
  if (a < b) { size_t tmp = a; a = b; b = tmp; }
  return a % 2 ? a * ((a - 1) / 2) + b : (a / 2) * (a - 1) + b;
}
static int valid_distance(double value) { return isnan(value) || (isfinite(value) && value >= 0); }

crpx_status crpx_distance_matrix_new(size_t size, crpx_distance_matrix_t *out)
{
  crpx_status status;
  size_t count;
  crpx_distance_matrix_t matrix;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  status = pair_count(size, &count);
  if (status != CRPX_OK) return status;
  if (count > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  matrix = calloc(1, sizeof(*matrix));
  if (!matrix) return CRPX_ERROR_ALLOCATION;
  matrix->values = count ? malloc(count * sizeof(double)) : NULL;
  if (count && !matrix->values) { free(matrix); return CRPX_ERROR_ALLOCATION; }
  matrix->size = size;
  matrix->count = count;
  for (size_t i = 0; i < count; ++i) matrix->values[i] = NAN;
  *out = matrix;
  return CRPX_OK;
}

void crpx_distance_matrix_free(crpx_distance_matrix_t matrix)
{
  if (matrix) { free(matrix->values); free(matrix); }
}
size_t crpx_distance_matrix_size(const struct crpx_distance_matrix_struct *matrix) { return matrix ? matrix->size : 0; }

crpx_status crpx_distance_matrix_clone(const struct crpx_distance_matrix_struct *source, crpx_distance_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!source) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_distance_matrix_new(source->size, out);
  if (status == CRPX_OK && source->count) memcpy((*out)->values, source->values, source->count * sizeof(double));
  return status;
}

crpx_status crpx_distance_matrix_set(crpx_distance_matrix_t matrix, size_t row, size_t column, double value)
{
  if (!matrix || !valid_distance(value)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (row >= matrix->size || column >= matrix->size) return CRPX_ERROR_OUT_OF_RANGE;
  if (row == column) return value == 0 ? CRPX_OK : CRPX_ERROR_INVALID_ARGUMENT;
  matrix->values[pair_index(row, column)] = value;
  return CRPX_OK;
}

crpx_status crpx_distance_matrix_get(const struct crpx_distance_matrix_struct *matrix, size_t row, size_t column,
                                   double *value)
{
  if (!matrix || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  if (row >= matrix->size || column >= matrix->size) return CRPX_ERROR_OUT_OF_RANGE;
  *value = row == column ? 0 : matrix->values[pair_index(row, column)];
  return CRPX_OK;
}

crpx_status crpx_distance_matrix_fill(crpx_distance_matrix_t matrix, double value)
{
  if (!matrix || !valid_distance(value)) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < matrix->count; ++i) matrix->values[i] = value;
  return CRPX_OK;
}

const double *crpx_distance_matrix_data(const struct crpx_distance_matrix_struct *matrix, size_t *count)
{
  if (count) *count = matrix ? matrix->count : 0;
  return matrix ? matrix->values : NULL;
}

crpx_status crpx_species_distances_new(size_t species, crpx_species_distances_t *out)
{
  crpx_status status;
  size_t count;
  crpx_species_distances_t distances;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  status = pair_count(species, &count);
  if (status != CRPX_OK) return status;
  if (count > SIZE_MAX / sizeof(struct species_pair)) return CRPX_ERROR_OVERFLOW;
  distances = calloc(1, sizeof(*distances));
  if (!distances) return CRPX_ERROR_ALLOCATION;
  distances->pairs = count ? calloc(count, sizeof(*distances->pairs)) : NULL;
  if (count && !distances->pairs) { free(distances); return CRPX_ERROR_ALLOCATION; }
  distances->size = species;
  distances->count = count;
  *out = distances;
  return CRPX_OK;
}

void crpx_species_distances_free(crpx_species_distances_t distances)
{
  if (distances) { free(distances->pairs); free(distances); }
}
void crpx_species_distances_clear(crpx_species_distances_t distances)
{
  if (distances && distances->count) memset(distances->pairs, 0, distances->count * sizeof(*distances->pairs));
}

crpx_status crpx_species_distances_add(crpx_species_distances_t distances, size_t a, size_t b, double value)
{
  if (!distances || !valid_distance(value)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (a >= distances->size || b >= distances->size) return CRPX_ERROR_OUT_OF_RANGE;
  if (a == b || isnan(value)) return CRPX_OK;
  struct species_pair *pair = distances->pairs + pair_index(a, b);
  if (pair->count == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  ++pair->count;
  pair->mean += (value - pair->mean) / (double)pair->count;
  if (pair->count == 1 || value < pair->minimum) pair->minimum = value;
  return CRPX_OK;
}

crpx_status crpx_species_distances_get(const struct crpx_species_distances_struct *distances, size_t a, size_t b,
                                     double *mean, double *minimum, uint64_t *count)
{
  if (!distances || !mean || !minimum || !count) return CRPX_ERROR_INVALID_ARGUMENT;
  if (a >= distances->size || b >= distances->size) return CRPX_ERROR_OUT_OF_RANGE;
  if (a == b) { *mean = *minimum = 0; *count = 0; return CRPX_OK; }
  const struct species_pair *pair = distances->pairs + pair_index(a, b);
  *mean = pair->count ? pair->mean : NAN;
  *minimum = pair->count ? pair->minimum : NAN;
  *count = pair->count;
  return CRPX_OK;
}

crpx_status crpx_species_distances_merge(crpx_species_distances_t to,
                                       const struct crpx_species_distances_struct *from)
{
  if (!to || !from || to->size != from->size) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < to->count; ++i)
    if (UINT64_MAX - to->pairs[i].count < from->pairs[i].count) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < to->count; ++i) {
    struct species_pair a = to->pairs[i], b = from->pairs[i];
    if (!b.count) continue;
    if (!a.count) { to->pairs[i] = b; continue; }
    to->pairs[i].count = a.count + b.count;
    to->pairs[i].mean = a.mean + (b.mean - a.mean) * ((double)b.count / (double)(a.count + b.count));
    to->pairs[i].minimum = a.minimum < b.minimum ? a.minimum : b.minimum;
  }
  return CRPX_OK;
}

crpx_status crpx_species_distances_add_genes(crpx_species_distances_t distances,
                                           const struct crpx_distance_matrix_struct *genes, const size_t *map)
{
  crpx_species_distances_t local;
  crpx_status status;
  if (!distances || !genes || (!map && genes->size)) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < genes->size; ++i)
    if (map[i] != SIZE_MAX && map[i] >= distances->size) return CRPX_ERROR_OUT_OF_RANGE;
  /* Accumulate separately so any error leaves the destination unchanged. */
  status = crpx_species_distances_new(distances->size, &local);
  if (status != CRPX_OK) return status;
  for (size_t i = 1; i < genes->size && status == CRPX_OK; ++i) if (map[i] != SIZE_MAX)
    for (size_t j = 0; j < i && status == CRPX_OK; ++j) if (map[j] != SIZE_MAX)
      status = crpx_species_distances_add(local, map[i], map[j], genes->values[pair_index(i, j)]);
  if (status == CRPX_OK) status = crpx_species_distances_merge(distances, local);
  crpx_species_distances_free(local);
  return status;
}

crpx_status crpx_species_distances_export(const struct crpx_species_distances_struct *distances,
                                        crpx_distance_matrix_t matrix, int use_minimum)
{
  if (!distances || !matrix || distances->size != matrix->size) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < matrix->count; ++i)
    matrix->values[i] = !distances->pairs[i].count ? NAN :
                       use_minimum ? distances->pairs[i].minimum : distances->pairs[i].mean;
  return CRPX_OK;
}
