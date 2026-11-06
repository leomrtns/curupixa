/* Reimplementation of biomcmc-lib lowlevel algorithms (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "lowlevel_algorithms.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

crpx_status crpx_assignment_solve(const double *costs, size_t rows, size_t columns, size_t *assignment, double *total)
{
  size_t *matched = NULL, *previous = NULL;
  long double *u = NULL, *v = NULL, *minimum = NULL;
  bool *used = NULL;
  crpx_status status = CRPX_OK;
  if (!total || rows > columns || (!costs && rows && columns) || (!assignment && rows))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (!rows) { *total = 0; return CRPX_OK; }
  if (columns == SIZE_MAX || rows == SIZE_MAX || columns > SIZE_MAX/rows ||
      columns+1 > SIZE_MAX/sizeof(long double) || columns+1 > SIZE_MAX/sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < rows*columns; ++i) if (!isfinite(costs[i]) || costs[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  u = calloc(rows+1, sizeof(*u)); v = calloc(columns+1, sizeof(*v));
  minimum = malloc((columns+1)*sizeof(*minimum));
  used = malloc((columns+1)*sizeof(*used));
  matched = calloc(columns+1, sizeof(*matched)); previous = malloc((columns+1)*sizeof(*previous));
  if (!u || !v || !minimum || !used || !matched || !previous) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t row = 1; row <= rows; ++row) {
    matched[0] = row;
    size_t column = 0;
    for (size_t j = 0; j <= columns; ++j) { minimum[j] = INFINITY; used[j] = false; }
    do {
      used[column] = true;
      size_t active_row = matched[column], next = 0;
      long double delta = INFINITY;
      for (size_t j = 1; j <= columns; ++j) if (!used[j]) {
        long double reduced = (long double)costs[(active_row-1)*columns+j-1]-u[active_row]-v[j];
        if (reduced < minimum[j]) { minimum[j] = reduced; previous[j] = column; }
        if (minimum[j] < delta) { delta = minimum[j]; next = j; }
      }
      if (!next || !isfinite(delta)) { status = CRPX_ERROR_OVERFLOW; goto done; }
      for (size_t j = 0; j <= columns; ++j) {
        if (used[j]) { u[matched[j]] += delta; v[j] -= delta; }
        else minimum[j] -= delta;
      }
      column = next;
    } while (matched[column]);
    do {
      size_t last = previous[column]; matched[column] = matched[last]; column = last;
    } while (column);
  }
  long double sum = 0;
  for (size_t j = 1; j <= columns; ++j) if (matched[j]) sum += costs[(matched[j]-1)*columns+j-1];
  if (!isfinite(sum) || sum > DBL_MAX) { status = CRPX_ERROR_OVERFLOW; goto done; }
  for (size_t j = 1; j <= columns; ++j) if (matched[j]) assignment[matched[j]-1] = j-1;
  *total = (double)sum;
done:
  free(u); free(v); free(minimum); free(used); free(matched); free(previous);
  return status;
}

struct edit_cost { uint64_t value; bool overflow; };
static struct edit_cost edit_add(struct edit_cost a, uint64_t b)
{
  if (a.overflow || a.value > UINT64_MAX-b) return (struct edit_cost){UINT64_MAX, true};
  return (struct edit_cost){a.value+b, false};
}
static struct edit_cost edit_min(struct edit_cost a, struct edit_cost b)
{
  return a.value < b.value || (a.value == b.value && !a.overflow) ? a : b;
}

crpx_status crpx_levenshtein_distance(const void *left, size_t n, const void *right, size_t m,
                                     uint64_t substitution, uint64_t indel, bool free_ends, uint64_t *distance)
{
  if ((!left && n) || (!right && m) || !distance) return CRPX_ERROR_INVALID_ARGUMENT;
  const unsigned char *a = left, *b = right;
  if (m > n) { const unsigned char *p = a; a = b; b = p; size_t size = n; n = m; m = size; }
  if (m == SIZE_MAX || m+1 > SIZE_MAX/sizeof(struct edit_cost)) return CRPX_ERROR_OVERFLOW;
  struct edit_cost *row = malloc((m+1)*sizeof(*row));
  if (!row) return CRPX_ERROR_ALLOCATION;
  row[0] = (struct edit_cost){0, false};
  for (size_t j = 1; j <= m; ++j) row[j] = edit_add(row[j-1], indel);
  struct edit_cost edge_minimum = row[m];
  for (size_t i = 1; i <= n; ++i) {
    struct edit_cost diagonal = row[0];
    row[0] = edit_add(row[0], indel);
    for (size_t j = 1; j <= m; ++j) {
      struct edit_cost before = row[j];
      struct edit_cost replace = a[i-1] == b[j-1] ? diagonal : edit_add(diagonal, substitution);
      struct edit_cost insert = edit_add(row[j-1], indel), remove = edit_add(before, indel);
      row[j] = edit_min(edit_min(replace, insert), remove);
      diagonal = before;
    }
    edge_minimum = edit_min(edge_minimum, row[m]);
  }
  struct edit_cost result = row[m];
  if (free_ends) {
    result = edge_minimum;
    for (size_t j = 0; j <= m; ++j) result = edit_min(result, row[j]);
  }
  free(row);
  if (result.overflow) return CRPX_ERROR_OVERFLOW;
  *distance = result.value;
  return CRPX_OK;
}

const char *crpx_string_last_substring(const char *text, const char *needle)
{
  if (!text || !needle) return NULL;
  if (!*needle) return text+strlen(text);
  const char *last = NULL, *found = text;
  while ((found = strstr(found, needle))) { last = found; ++found; }
  return last;
}
size_t crpx_string_common_prefix(const char *a, const char *b)
{
  size_t count = 0;
  if (a && b) while (a[count] && a[count] == b[count]) ++count;
  return count;
}
