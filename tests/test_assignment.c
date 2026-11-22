/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
#include "assignment.h"
#include "lowlevel_algorithms.h"
#include "global/random_number.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); return 1; } } while (0)

static double
brute (const double *costs, size_t rows, size_t columns, size_t row, unsigned used)
{
  if (row == rows) return 0;
  double best = INFINITY;
  for (size_t j = 0; j < columns; j++) if (!(used & (1U << j))) {
    double candidate = costs[row * columns + j] + brute (costs, rows, columns, row + 1, used | (1U << j));
    if (candidate < best) best = candidate;
  }
  return best;
}

int
main (void)
{
  crpx_assignment_workspace *work = NULL;
  CHECK (crpx_assignment_workspace_new (4, 6, &work) == CRPX_OK);
  crpx_rng_stream stream;
  crpx_rng_stream_init (&stream, 378, 0);
  double costs[24], copy[24], total, old_total;
  size_t assignment[4], old_assignment[4];
  for (size_t rows = 1; rows <= 4; rows++) for (size_t columns = rows; columns <= 6; columns++)
    for (size_t repeat = 0; repeat < 40; repeat++) {
      for (size_t i = 0; i < rows * columns; i++) costs[i] = (double)(crpx_rng_stream_next (&stream) % 13);
      memcpy (copy, costs, sizeof (costs));
      CHECK (crpx_assignment_jv_with_workspace (work, costs, rows, columns, assignment, &total) == CRPX_OK);
      CHECK (crpx_assignment_solve (costs, rows, columns, old_assignment, &old_total) == CRPX_OK);
      CHECK (total == old_total && total == brute (costs, rows, columns, 0, 0));
      CHECK (!memcmp (copy, costs, sizeof (costs)));
      unsigned used = 0;
      double actual = 0;
      for (size_t i = 0; i < rows; i++) {
        CHECK (assignment[i] < columns && !(used & (1U << assignment[i])));
        used |= 1U << assignment[i]; actual += costs[i * columns + assignment[i]];
      }
      CHECK (actual == total);
    }
  memset (costs, 0, sizeof (costs));
  CHECK (crpx_assignment_jv (costs, 4, 6, assignment, &total) == CRPX_OK && total == 0);
  for (size_t i = 0; i < 4; i++) CHECK (assignment[i] == i);
  total = -1; assignment[0] = 99;
  costs[0] = NAN;
  CHECK (crpx_assignment_jv_with_workspace (work, costs, 2, 2, assignment, &total) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK (assignment[0] == 99 && total == -1);
  for (size_t i = 0; i < 4; i++) costs[i] = DBL_MAX;
  CHECK (crpx_assignment_jv (costs, 2, 2, assignment, &total) == CRPX_ERROR_OVERFLOW);
  CHECK (assignment[0] == 99 && total == -1);
  costs[0] = costs[3] = 0;
  CHECK (crpx_assignment_jv (costs, 2, 2, assignment, &total) == CRPX_OK && total == 0);
  CHECK (crpx_assignment_jv_with_workspace (work, costs, 5, 5, assignment, &total) == CRPX_ERROR_OUT_OF_RANGE);
  CHECK (crpx_assignment_jv (NULL, 0, 0, NULL, &total) == CRPX_OK && total == 0);
  CHECK (crpx_assignment_jv (costs, 3, 2, assignment, &total) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_assignment_workspace_free (work);
  return 0;
}
