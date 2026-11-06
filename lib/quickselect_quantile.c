/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "quickselect_quantile.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void swap_value(double *a, double *b)
{
  double value = *a;
  *a = *b;
  *b = value;
}

static double select_value(double *values, size_t size, size_t rank)
{
  size_t left = 0, right = size;
  while (right - left > 1) {
    double pivot = values[left + (right - left) / 2];
    size_t lower = left, current = left, upper = right;
    while (current < upper) {
      if (values[current] < pivot) swap_value(values + lower++, values + current++);
      else if (values[current] > pivot) swap_value(values + current, values + --upper);
      else current++;
    }
    if (rank < lower) right = lower;
    else if (rank >= upper) left = upper;
    else return values[rank];
  }
  return values[rank];
}

crpx_status crpx_select_double(double *values, size_t size, size_t rank, double *out)
{
  if (!values || !size || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (rank >= size) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t i = 0; i < size; ++i) if (isnan(values[i])) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = select_value(values, size, rank);
  return CRPX_OK;
}

crpx_status crpx_quantile_vector_double(const double *values, size_t size, const double *probabilities,
                                      size_t count, double *out)
{
  if (!values || !size || (count && (!probabilities || !out))) return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < size; ++i) if (isnan(values[i])) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < count; ++i)
    if (!isfinite(probabilities[i]) || probabilities[i] < 0 || probabilities[i] > 1)
      return CRPX_ERROR_INVALID_ARGUMENT;
  if (!count) return CRPX_OK;
  double *copy = malloc(size * sizeof(double));
  if (!copy) return CRPX_ERROR_ALLOCATION;
  memcpy(copy, values, size * sizeof(double));
  for (size_t i = 0; i < count; ++i) {
    long double position = (long double)size * probabilities[i];
    size_t rank = position >= (long double)size ? size - 1 : (size_t)position;
    out[i] = select_value(copy, size, rank);
  }
  free(copy);
  return CRPX_OK;
}

crpx_status crpx_quantile_double(const double *values, size_t size, double probability, double *out)
{
  return crpx_quantile_vector_double(values, size, &probability, 1, out);
}
