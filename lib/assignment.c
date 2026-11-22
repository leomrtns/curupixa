/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT AND BSD-3-Clause */
/*! \file assignment.c
    Modified Jonker-Volgenant shortest augmenting paths, following Crouse (2016), doi:10.1109/TAES.2016.140952.
    C17 adaptation of PM Larsen's BSD-licensed SciPy rectangular_lsap.cpp; SciPy notices are retained in
    LICENSES/Assignment-SciPy-BSD-3-Clause.txt. Changes: reusable arrays, long-double intermediates, checked C API. */
#include "assignment.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t augmenting_path (crpx_assignment_workspace *work, const double *costs, size_t rows,
                               size_t columns, size_t row, long double *minimum);

crpx_status
crpx_assignment_workspace_new (size_t rows, size_t columns, crpx_assignment_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (rows > columns) return CRPX_ERROR_INVALID_ARGUMENT;
  if (columns > SIZE_MAX / sizeof (long double) || columns > SIZE_MAX / sizeof (size_t)) return CRPX_ERROR_OVERFLOW;
  crpx_assignment_workspace *work = calloc (1, sizeof (*work));
  if (!work) return CRPX_ERROR_ALLOCATION;
  work->rows = rows; work->columns = columns;
  if (rows) {
    work->row_potential = malloc (rows * sizeof (*work->row_potential));
    work->row_column = malloc (rows * sizeof (*work->row_column));
    work->visited_rows = malloc (rows * sizeof (*work->visited_rows));
  }
  if (columns) {
    work->column_potential = malloc (columns * sizeof (*work->column_potential));
    work->distance = malloc (columns * sizeof (*work->distance));
    work->column_row = malloc (columns * sizeof (*work->column_row));
    work->previous = malloc (columns * sizeof (*work->previous));
    work->remaining = malloc (columns * sizeof (*work->remaining));
    work->visited_columns = malloc (columns * sizeof (*work->visited_columns));
  }
  if ((rows && (!work->row_potential || !work->row_column || !work->visited_rows)) ||
      (columns && (!work->column_potential || !work->distance || !work->column_row || !work->previous ||
                   !work->remaining || !work->visited_columns))) {
    crpx_assignment_workspace_free (work); return CRPX_ERROR_ALLOCATION;
  }
  *out = work;
  return CRPX_OK;
}

void
crpx_assignment_workspace_free (crpx_assignment_workspace *work)
{
  if (!work) return;
  free (work->row_potential); free (work->column_potential); free (work->distance);
  free (work->row_column); free (work->column_row); free (work->previous); free (work->remaining);
  free (work->visited_rows); free (work->visited_columns); free (work);
}

crpx_status
crpx_assignment_jv (const double *costs, size_t rows, size_t columns, size_t *assignment, double *total)
{ // One-off convenience call; retain a workspace for repeated tree comparisons or matching problems.
  if (!total || rows > columns || (!costs && rows) || (!assignment && rows)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!rows) { *total = 0; return CRPX_OK; }
  if (columns > SIZE_MAX / rows) return CRPX_ERROR_OVERFLOW;
  crpx_assignment_workspace *work = NULL;
  crpx_status status = crpx_assignment_workspace_new (rows, columns, &work);
  if (status == CRPX_OK) status = crpx_assignment_jv_with_workspace (work, costs, rows, columns, assignment, total);
  crpx_assignment_workspace_free (work);
  return status;
}

crpx_status
crpx_assignment_jv_with_workspace (crpx_assignment_workspace *work, const double *costs,
                                   size_t rows, size_t columns, size_t *assignment, double *total)
{ // Reinitialise duals and matching, then augment one row at a time without allocating.
  if (!work || !total || rows > columns || (!costs && rows) || (!assignment && rows)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (rows > work->rows || columns > work->columns) return CRPX_ERROR_OUT_OF_RANGE;
  if (!rows) { *total = 0; return CRPX_OK; }
  if (columns > SIZE_MAX / rows) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < rows * columns; i++)
    if (!isfinite (costs[i]) || costs[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < rows; i++) { work->row_potential[i] = 0; work->row_column[i] = SIZE_MAX; }
  for (size_t j = 0; j < columns; j++) { work->column_potential[j] = 0; work->column_row[j] = SIZE_MAX; }
  for (size_t row = 0; row < rows; row++) {
    long double minimum;
    size_t sink = augmenting_path (work, costs, rows, columns, row, &minimum);
    if (sink == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
    work->row_potential[row] += minimum;
    for (size_t i = 0; i < rows; i++) if (work->visited_rows[i] && i != row)
      work->row_potential[i] += minimum - work->distance[work->row_column[i]];
    for (size_t j = 0; j < columns; j++) if (work->visited_columns[j])
      work->column_potential[j] -= minimum - work->distance[j];
    size_t column = sink;
    for (;;) {
      size_t current = work->previous[column], previous = work->row_column[current];
      work->column_row[column] = current; work->row_column[current] = column;
      if (current == row) break;
      column = previous;
    }
  }
  long double sum = 0;
  for (size_t i = 0; i < rows; i++) sum += costs[i * columns + work->row_column[i]];
  if (!isfinite (sum) || sum > DBL_MAX) return CRPX_ERROR_OVERFLOW;
  memcpy (assignment, work->row_column, rows * sizeof (*assignment));
  *total = (double)sum;
  return CRPX_OK;
}

static size_t
augmenting_path (crpx_assignment_workspace *work, const double *costs, size_t rows,
                  size_t columns, size_t row, long double *minimum)
{
  memset (work->visited_rows, 0, rows);
  memset (work->visited_columns, 0, columns);
  for (size_t j = 0; j < columns; j++) { work->remaining[j] = columns - j - 1; work->distance[j] = INFINITY; }
  size_t remaining = columns;
  long double current_minimum = 0;
  for (;;) {
    size_t index = SIZE_MAX;
    long double lowest = INFINITY;
    work->visited_rows[row] = 1;
    for (size_t i = 0; i < remaining; i++) {
      size_t column = work->remaining[i];
      long double reduced = current_minimum + costs[row * columns + column] -
                            work->row_potential[row] - work->column_potential[column];
      if (reduced < work->distance[column]) { work->distance[column] = reduced; work->previous[column] = row; }
      if (work->distance[column] < lowest ||
          (work->distance[column] == lowest && work->column_row[column] == SIZE_MAX)) {
        lowest = work->distance[column]; index = i;
      }
    }
    if (index == SIZE_MAX || !isfinite (lowest)) return SIZE_MAX;
    current_minimum = lowest;
    size_t column = work->remaining[index];
    work->visited_columns[column] = 1;
    work->remaining[index] = work->remaining[--remaining];
    if (work->column_row[column] == SIZE_MAX) { *minimum = current_minimum; return column; }
    row = work->column_row[column];
  }
}
