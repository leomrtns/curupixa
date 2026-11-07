/* Reimplementation of biomcmc-lib distance_matrix (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_DISTANCE_MATRIX_H
#define CURUPIXA_DISTANCE_MATRIX_H
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_distance_matrix_struct *crpx_distance_matrix_t;
typedef struct crpx_species_distances_struct *crpx_species_distances_t;

/* Packed symmetric off-diagonal storage, n(n-1)/2 doubles. Diagonal is zero.
 * Missing values are NaN; inputs must be non-negative finite values or NaN. */
crpx_status crpx_distance_matrix_new(size_t size, crpx_distance_matrix_t *out);
void crpx_distance_matrix_free(crpx_distance_matrix_t matrix);
size_t crpx_distance_matrix_size(const struct crpx_distance_matrix_struct *matrix);
crpx_status crpx_distance_matrix_clone(const struct crpx_distance_matrix_struct *source, crpx_distance_matrix_t *out);
crpx_status crpx_distance_matrix_set(crpx_distance_matrix_t matrix, size_t row, size_t column, double value);
crpx_status crpx_distance_matrix_get(const struct crpx_distance_matrix_struct *matrix, size_t row, size_t column,
                                   double *value);
crpx_status crpx_distance_matrix_fill(crpx_distance_matrix_t matrix, double value);
const double *crpx_distance_matrix_data(const struct crpx_distance_matrix_struct *matrix, size_t *count);

/* Keep means, minima and counts separately; no overloaded matrix triangles.
 * Each observed pair contributes once. NaN observations are skipped. */
crpx_status crpx_species_distances_new(size_t species, crpx_species_distances_t *out);
void crpx_species_distances_free(crpx_species_distances_t distances);
void crpx_species_distances_clear(crpx_species_distances_t distances);
crpx_status crpx_species_distances_add(crpx_species_distances_t distances, size_t a, size_t b, double value);
crpx_status crpx_species_distances_get(const struct crpx_species_distances_struct *distances, size_t a, size_t b,
                                     double *mean, double *minimum, uint64_t *count);
crpx_status crpx_species_distances_merge(crpx_species_distances_t to,
                                       const struct crpx_species_distances_struct *from);
/* map[i] gives the species for gene i; SIZE_MAX skips a gene. Validates the full
 * map before adding any observation. */
crpx_status crpx_species_distances_add_genes(crpx_species_distances_t distances,
                                           const struct crpx_distance_matrix_struct *genes, const size_t *map);
crpx_status crpx_species_distances_export(const struct crpx_species_distances_struct *distances,
                                        crpx_distance_matrix_t matrix, int use_minimum);

#ifdef __cplusplus
}
#endif
#endif
