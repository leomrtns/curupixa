/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file assignment.h
    Rectangular Jonker-Volgenant assignment alongside the Hungarian solver in lowlevel_algorithms.h.
    The workspace avoids allocations when successive cost matrices fit the original dimensions. */
#ifndef CURUPIXA_ASSIGNMENT_H
#define CURUPIXA_ASSIGNMENT_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Construct with workspace_new; use each workspace in one calculation at a time. Dimensions are capacities. */
typedef struct {
  size_t rows, columns;
  long double *row_potential, *column_potential, *distance;
  size_t *row_column, *column_row, *previous, *remaining;
  unsigned char *visited_rows, *visited_columns;
} crpx_assignment_workspace;

crpx_status crpx_assignment_workspace_new (size_t rows, size_t columns, crpx_assignment_workspace **out);
void crpx_assignment_workspace_free (crpx_assignment_workspace *workspace);
/* Same matrix and output conventions as crpx_assignment_solve: rows<=columns, finite non-negative row-major costs.
   Outputs remain unchanged on failure. Equal-cost assignments can differ between the two algorithms.
   Time O(rows^2*columns), scratch O(rows+columns); the input matrix is unchanged. */
crpx_status crpx_assignment_jv (const double *costs, size_t rows, size_t columns, size_t *assignment, double *total);
crpx_status crpx_assignment_jv_with_workspace (crpx_assignment_workspace *workspace, const double *costs,
                                              size_t rows, size_t columns, size_t *assignment, double *total);
#ifdef __cplusplus
}
#endif
#endif
