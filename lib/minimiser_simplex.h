/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_MINIMISER_SIMPLEX_H
#define CURUPIXA_MINIMISER_SIMPLEX_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_simplex crpx_simplex;
typedef crpx_status (*crpx_objective)(const double *, size_t, void *, double *);
/* Bounded Nelder-Mead; one object per worker. Callback errors/non-finite values
 * return to the caller without committing a partial iteration. Callback receives
 * point, dimensions, caller-owned parameters and output value. */
crpx_status crpx_simplex_new(size_t dimensions, crpx_objective objective, void *parameters, crpx_simplex **out);
void crpx_simplex_free(crpx_simplex *simplex);
crpx_status crpx_simplex_initialise(crpx_simplex *simplex, const double *point, const double *steps);
crpx_status crpx_simplex_iterate(crpx_simplex *simplex);
/* Both absolute tolerances must hold: maximum vertex distance from the best
 * point and maximum objective difference. CONVERGENCE means the budget expired;
 * the best point/value remain available. */
crpx_status crpx_simplex_run(crpx_simplex *simplex, size_t maximum_iterations, double point_tolerance,
                            double value_tolerance, size_t *iterations);
const double *crpx_simplex_best_point(const crpx_simplex *simplex);
double crpx_simplex_best_value(const crpx_simplex *simplex);
double crpx_simplex_size(const crpx_simplex *simplex);
size_t crpx_simplex_dimensions(const crpx_simplex *simplex);
#ifdef __cplusplus
}
#endif
#endif
