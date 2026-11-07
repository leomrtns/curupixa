/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_SVD_H
#define CURUPIXA_SVD_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_matrix_struct *crpx_matrix_t;
typedef struct crpx_sparse_matrix_struct *crpx_sparse_matrix_t;
typedef struct crpx_svd_struct *crpx_svd_t;

/* Matrices own their storage. Dense arrays use row-major order. Concurrent reads
 * are safe; callers serialise changes and destruction of the same matrix. */
crpx_status crpx_matrix_new(size_t rows, size_t columns, crpx_matrix_t *out);
crpx_status crpx_matrix_from_array(const double *values, size_t rows, size_t columns, crpx_matrix_t *out);
void crpx_matrix_free(crpx_matrix_t matrix);
size_t crpx_matrix_rows(const struct crpx_matrix_struct *matrix);
size_t crpx_matrix_columns(const struct crpx_matrix_struct *matrix);
double *crpx_matrix_data(crpx_matrix_t matrix);
const double *crpx_matrix_const_data(const struct crpx_matrix_struct *matrix);
crpx_status crpx_matrix_get(const struct crpx_matrix_struct *matrix, size_t row, size_t column, double *out);
crpx_status crpx_matrix_set(crpx_matrix_t matrix, size_t row, size_t column, double value);
crpx_status crpx_matrix_transpose(const struct crpx_matrix_struct *matrix, crpx_matrix_t *out);
/* Input and output vectors must not overlap. Lengths are checked explicitly. */
crpx_status crpx_matrix_multiply_vector(const struct crpx_matrix_struct *matrix, const double *input,
                                      size_t input_size, double *output, size_t output_size);

/* Compressed sparse columns: offsets has columns+1 entries, starts at zero,
 * ends at nonzero, and is non-decreasing. Row indices must be strictly increasing
 * within each column. Construction copies all arrays. */
crpx_status crpx_sparse_matrix_new(size_t rows, size_t columns, size_t nonzero, const size_t *offsets,
                                  const size_t *row_indices, const double *values, crpx_sparse_matrix_t *out);
void crpx_sparse_matrix_free(crpx_sparse_matrix_t matrix);
size_t crpx_sparse_matrix_rows(const struct crpx_sparse_matrix_struct *matrix);
size_t crpx_sparse_matrix_columns(const struct crpx_sparse_matrix_struct *matrix);
size_t crpx_sparse_matrix_nonzero(const struct crpx_sparse_matrix_struct *matrix);
const size_t *crpx_sparse_matrix_offsets(const struct crpx_sparse_matrix_struct *matrix);
const size_t *crpx_sparse_matrix_row_indices(const struct crpx_sparse_matrix_struct *matrix);
const double *crpx_sparse_matrix_values(const struct crpx_sparse_matrix_struct *matrix);
crpx_status crpx_sparse_matrix_to_dense(const struct crpx_sparse_matrix_struct *matrix, crpx_matrix_t *out);
crpx_status crpx_sparse_matrix_from_dense(const struct crpx_matrix_struct *matrix, double zero_tolerance,
                                         crpx_sparse_matrix_t *out);
crpx_status crpx_sparse_matrix_transpose(const struct crpx_sparse_matrix_struct *matrix, crpx_sparse_matrix_t *out);
crpx_status crpx_sparse_matrix_multiply_vector(const struct crpx_sparse_matrix_struct *matrix, const double *input,
                                             size_t input_size, double *output, size_t output_size);

/* Thin SVD A = U diag(s) Vt. Singular values decrease. U has rows*rank entries,
 * Vt has rank*columns entries. rank=0 keeps min(rows,columns). tolerance=0 uses
 * 8*DBL_EPSILON; sweeps=0 uses 100. All workspaces are local and failures leave
 * *out=NULL. CRPX_ERROR_CONVERGENCE reports an exhausted sweep limit.
 * This dependency-free dense Jacobi backend targets small/moderate matrices;
 * sparse input is explicitly materialised, so storage is O(rows*columns).
 * Columns for exact zero singular values may be zero vectors. */
crpx_status crpx_svd_decompose(const struct crpx_matrix_struct *matrix, size_t rank, double tolerance,
                              size_t sweeps, crpx_svd_t *out);
crpx_status crpx_svd_decompose_sparse(const struct crpx_sparse_matrix_struct *matrix, size_t rank,
                                     double tolerance, size_t sweeps, crpx_svd_t *out);
void crpx_svd_free(crpx_svd_t decomposition);
size_t crpx_svd_rank(const struct crpx_svd_struct *decomposition);
size_t crpx_svd_rows(const struct crpx_svd_struct *decomposition);
size_t crpx_svd_columns(const struct crpx_svd_struct *decomposition);
const double *crpx_svd_singular_values(const struct crpx_svd_struct *decomposition);
const double *crpx_svd_left_vectors(const struct crpx_svd_struct *decomposition);
const double *crpx_svd_right_vectors_transposed(const struct crpx_svd_struct *decomposition);
crpx_status crpx_svd_reconstruct(const struct crpx_svd_struct *decomposition, crpx_matrix_t *out);

#ifdef __cplusplus
}
#endif
#endif
