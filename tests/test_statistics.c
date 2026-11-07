/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "quickselect_quantile.h"
#include "empirical_frequency.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; } } while (0)

static int compare_double(const void *a, const void *b)
{ return (*(const double *)a > *(const double *)b) - (*(const double *)a < *(const double *)b); }

int main(void)
{
  for (size_t size = 1; size < 100; ++size) {
    double original[100], sorted[100], copy[100], result;
    for (size_t i = 0; i < size; ++i) original[i] = (double)((i * 37 + size * 13) % 17);
    memcpy(sorted, original, size * sizeof(double));
    qsort(sorted, size, sizeof(double), compare_double);
    for (size_t rank = 0; rank < size; ++rank) {
      memcpy(copy, original, size * sizeof(double));
      CHECK(crpx_select_double(copy, size, rank, &result) == CRPX_OK);
      CHECK(result == sorted[rank]);
    }
    CHECK(crpx_quantile_double(original, size, 1, &result) == CRPX_OK);
    CHECK(result == sorted[size - 1]);
    CHECK(crpx_quantile_double(original, size, 0, &result) == CRPX_OK);
    CHECK(result == sorted[0]);
    CHECK(crpx_quantile_double(original, size, 0.5, &result) == CRPX_OK);
    CHECK(result == sorted[size / 2]);
  }
  double invalid[] = {1, NAN}, result = 42;
  CHECK(crpx_select_double(invalid, 2, 0, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(result == 42 && invalid[0] == 1 && isnan(invalid[1]));
  CHECK(crpx_quantile_double(invalid, 2, NAN, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(crpx_quantile_double(invalid, 0, 0, &result) == CRPX_ERROR_INVALID_ARGUMENT);

  int64_t values[] = {INT64_MIN, 5, INT64_MAX, 5, -1, -1};
  uint64_t weights[] = {0, 2, 1, 3, 2, 3};
  crpx_empfreq_t frequencies = NULL, merged = NULL, empty = NULL;
  CHECK(crpx_empfreq_new(values, weights, 6, &frequencies) == CRPX_OK);
  CHECK(crpx_empfreq_size(frequencies) == 3 && crpx_empfreq_total(frequencies) == 11);
  const crpx_empfreq_entry *entries = crpx_empfreq_entries(frequencies);
  CHECK(entries[0].value == -1 && entries[0].frequency == 5);
  CHECK(entries[1].value == 5 && entries[1].frequency == 5);
  int64_t mode, minimum, maximum;
  CHECK(crpx_empfreq_mode(frequencies, &mode) == CRPX_OK && mode == -1);
  CHECK(crpx_empfreq_range(frequencies, &minimum, &maximum) == CRPX_OK);
  CHECK(minimum == -1 && maximum == INT64_MAX);
  CHECK(crpx_empfreq_sort(frequencies, false) == CRPX_OK);
  CHECK(crpx_empfreq_entries(frequencies)[0].value == INT64_MAX);
  CHECK(crpx_empfreq_mode(frequencies, &mode) == CRPX_OK && mode == -1);
  CHECK(crpx_empfreq_merge(frequencies, frequencies, &merged) == CRPX_OK);
  CHECK(crpx_empfreq_total(merged) == 22 && crpx_empfreq_entries(merged)[0].frequency == 10);
  crpx_empfreq_free(merged);
  CHECK(crpx_empfreq_new(NULL, NULL, 0, &empty) == CRPX_OK);
  CHECK(crpx_empfreq_mode(empty, &mode) == CRPX_ERROR_NOT_FOUND);
  CHECK(crpx_empfreq_merge(frequencies, empty, &merged) == CRPX_OK);
  CHECK(crpx_empfreq_total(merged) == 11);
  crpx_empfreq_free(merged);
  uint64_t huge_weights[] = {UINT64_MAX, 1};
  CHECK(crpx_empfreq_new(values, huge_weights, 2, &merged) == CRPX_ERROR_OVERFLOW && !merged);
  crpx_empfreq_free(empty);
  crpx_empfreq_free(frequencies);

  size_t indices[6];
  CHECK(crpx_rank_int64(values, 6, false, indices) == CRPX_OK);
  CHECK(indices[0] == 0 && indices[1] == 4 && indices[2] == 5 && indices[5] == 2);
  CHECK(crpx_rank_int64(values, 6, true, indices) == CRPX_OK);
  CHECK(indices[0] == 2 && indices[1] == 1 && indices[2] == 3 && indices[5] == 0);
  double doubles[] = {INFINITY, -INFINITY, 0, -0.0};
  CHECK(crpx_rank_double(doubles, 4, false, indices) == CRPX_OK);
  CHECK(indices[0] == 1 && indices[1] == 2 && indices[2] == 3 && indices[3] == 0);
  CHECK(crpx_rank_double(invalid, 2, false, indices) == CRPX_ERROR_INVALID_ARGUMENT);
  return 0;
}
