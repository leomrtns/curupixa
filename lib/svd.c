/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "svd.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct crpx_matrix_struct { size_t rows, columns; double *values; };
struct crpx_sparse_matrix_struct {
  size_t rows, columns, nonzero;
  size_t *offsets, *row_indices;
  double *values;
};
struct crpx_svd_struct { size_t rows, columns, rank; double *singular, *left, *right; };

static bool product_fits(size_t a, size_t b, size_t item_size)
{ return !a || b <= SIZE_MAX / a / item_size; }

crpx_status crpx_matrix_new(size_t rows, size_t columns, crpx_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!rows || !columns) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!product_fits(rows, columns, sizeof(double))) return CRPX_ERROR_OVERFLOW;
  crpx_matrix_t matrix = calloc(1, sizeof(*matrix));
  if (!matrix) return CRPX_ERROR_ALLOCATION;
  matrix->values = calloc(rows * columns, sizeof(double));
  if (!matrix->values) { free(matrix); return CRPX_ERROR_ALLOCATION; }
  matrix->rows = rows;
  matrix->columns = columns;
  *out = matrix;
  return CRPX_OK;
}

crpx_status crpx_matrix_from_array(const double *values, size_t rows, size_t columns, crpx_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!values) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_matrix_new(rows, columns, out);
  if (status == CRPX_OK) memcpy((*out)->values, values, rows * columns * sizeof(double));
  return status;
}

void crpx_matrix_free(crpx_matrix_t matrix) { if (matrix) { free(matrix->values); free(matrix); } }
size_t crpx_matrix_rows(const struct crpx_matrix_struct *matrix) { return matrix ? matrix->rows : 0; }
size_t crpx_matrix_columns(const struct crpx_matrix_struct *matrix) { return matrix ? matrix->columns : 0; }
double *crpx_matrix_data(crpx_matrix_t matrix) { return matrix ? matrix->values : NULL; }
const double *crpx_matrix_const_data(const struct crpx_matrix_struct *matrix) { return matrix ? matrix->values : NULL; }

crpx_status crpx_matrix_get(const struct crpx_matrix_struct *matrix, size_t row, size_t column, double *out)
{
  if (!matrix || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (row >= matrix->rows || column >= matrix->columns) return CRPX_ERROR_OUT_OF_RANGE;
  *out = matrix->values[row * matrix->columns + column];
  return CRPX_OK;
}

crpx_status crpx_matrix_set(crpx_matrix_t matrix, size_t row, size_t column, double value)
{
  if (!matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  if (row >= matrix->rows || column >= matrix->columns) return CRPX_ERROR_OUT_OF_RANGE;
  matrix->values[row * matrix->columns + column] = value;
  return CRPX_OK;
}

crpx_status crpx_matrix_transpose(const struct crpx_matrix_struct *matrix, crpx_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_matrix_new(matrix->columns, matrix->rows, out);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < matrix->rows; ++i)
    for (size_t j = 0; j < matrix->columns; ++j)
      (*out)->values[j * matrix->rows + i] = matrix->values[i * matrix->columns + j];
  return CRPX_OK;
}

crpx_status crpx_matrix_multiply_vector(const struct crpx_matrix_struct *matrix, const double *input,
                                      size_t input_size, double *output, size_t output_size)
{
  if (!matrix || !input || !output) return CRPX_ERROR_INVALID_ARGUMENT;
  if (input_size != matrix->columns || output_size != matrix->rows) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < matrix->rows; ++i) {
    long double sum = 0;
    for (size_t j = 0; j < matrix->columns; ++j) sum += (long double)matrix->values[i * matrix->columns + j] * input[j];
    output[i] = (double)sum;
  }
  return CRPX_OK;
}

crpx_status crpx_sparse_matrix_new(size_t rows, size_t columns, size_t nonzero, const size_t *offsets,
                                  const size_t *row_indices, const double *values, crpx_sparse_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!rows || !columns || !offsets || (nonzero && (!row_indices || !values))) return CRPX_ERROR_INVALID_ARGUMENT;
  if (columns == SIZE_MAX || columns + 1 > SIZE_MAX / sizeof(size_t) ||
      nonzero > SIZE_MAX / sizeof(size_t) || nonzero > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  if (offsets[0] || offsets[columns] != nonzero) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t j = 0; j < columns; ++j) {
    if (offsets[j] > offsets[j + 1] || offsets[j + 1] > nonzero) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t k = offsets[j]; k < offsets[j + 1]; ++k)
      if (row_indices[k] >= rows || (k > offsets[j] && row_indices[k] <= row_indices[k - 1]))
        return CRPX_ERROR_INVALID_ARGUMENT;
  }
  crpx_sparse_matrix_t matrix = calloc(1, sizeof(*matrix));
  if (!matrix) return CRPX_ERROR_ALLOCATION;
  matrix->offsets = malloc((columns + 1) * sizeof(size_t));
  if (nonzero) {
    matrix->row_indices = malloc(nonzero * sizeof(size_t));
    matrix->values = malloc(nonzero * sizeof(double));
  }
  if (!matrix->offsets || (nonzero && (!matrix->row_indices || !matrix->values))) {
    crpx_sparse_matrix_free(matrix);
    return CRPX_ERROR_ALLOCATION;
  }
  matrix->rows = rows; matrix->columns = columns; matrix->nonzero = nonzero;
  memcpy(matrix->offsets, offsets, (columns + 1) * sizeof(size_t));
  if (nonzero) {
    memcpy(matrix->row_indices, row_indices, nonzero * sizeof(size_t));
    memcpy(matrix->values, values, nonzero * sizeof(double));
  }
  *out = matrix;
  return CRPX_OK;
}

void crpx_sparse_matrix_free(crpx_sparse_matrix_t matrix)
{
  if (matrix) { free(matrix->offsets); free(matrix->row_indices); free(matrix->values); free(matrix); }
}
size_t crpx_sparse_matrix_rows(const struct crpx_sparse_matrix_struct *matrix) { return matrix ? matrix->rows : 0; }
size_t crpx_sparse_matrix_columns(const struct crpx_sparse_matrix_struct *matrix) { return matrix ? matrix->columns : 0; }
size_t crpx_sparse_matrix_nonzero(const struct crpx_sparse_matrix_struct *matrix) { return matrix ? matrix->nonzero : 0; }
const size_t *crpx_sparse_matrix_offsets(const struct crpx_sparse_matrix_struct *matrix)
{ return matrix ? matrix->offsets : NULL; }
const size_t *crpx_sparse_matrix_row_indices(const struct crpx_sparse_matrix_struct *matrix)
{ return matrix ? matrix->row_indices : NULL; }
const double *crpx_sparse_matrix_values(const struct crpx_sparse_matrix_struct *matrix)
{ return matrix ? matrix->values : NULL; }

crpx_status crpx_sparse_matrix_to_dense(const struct crpx_sparse_matrix_struct *matrix, crpx_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_matrix_new(matrix->rows, matrix->columns, out);
  if (status != CRPX_OK) return status;
  for (size_t j = 0; j < matrix->columns; ++j)
    for (size_t k = matrix->offsets[j]; k < matrix->offsets[j + 1]; ++k)
      (*out)->values[matrix->row_indices[k] * matrix->columns + j] = matrix->values[k];
  return CRPX_OK;
}

crpx_status crpx_sparse_matrix_from_dense(const struct crpx_matrix_struct *matrix, double zero_tolerance,
                                         crpx_sparse_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix || !isfinite(zero_tolerance) || zero_tolerance < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  if (matrix->columns == SIZE_MAX || matrix->columns + 1 > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  size_t nonzero = 0;
  for (size_t i = 0; i < matrix->rows * matrix->columns; ++i) {
    if (!isfinite(matrix->values[i])) return CRPX_ERROR_INVALID_ARGUMENT;
    nonzero += fabs(matrix->values[i]) > zero_tolerance;
  }
  if (nonzero > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  size_t *offsets = calloc(matrix->columns + 1, sizeof(size_t));
  size_t *indices = nonzero ? malloc(nonzero * sizeof(size_t)) : NULL;
  double *values = nonzero ? malloc(nonzero * sizeof(double)) : NULL;
  if (!offsets || (nonzero && (!indices || !values))) {
    free(offsets); free(indices); free(values); return CRPX_ERROR_ALLOCATION;
  }
  size_t count = 0;
  for (size_t j = 0; j < matrix->columns; ++j) {
    for (size_t i = 0; i < matrix->rows; ++i)
      if (fabs(matrix->values[i * matrix->columns + j]) > zero_tolerance) {
        indices[count] = i;
        values[count++] = matrix->values[i * matrix->columns + j];
      }
    offsets[j + 1] = count;
  }
  crpx_status status = crpx_sparse_matrix_new(matrix->rows, matrix->columns, nonzero, offsets, indices, values, out);
  free(offsets); free(indices); free(values);
  return status;
}

crpx_status crpx_sparse_matrix_transpose(const struct crpx_sparse_matrix_struct *matrix, crpx_sparse_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix) return CRPX_ERROR_INVALID_ARGUMENT;
  if (matrix->rows == SIZE_MAX || matrix->rows + 1 > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  size_t *offsets = calloc(matrix->rows + 1, sizeof(size_t));
  size_t *next = calloc(matrix->rows, sizeof(size_t));
  size_t *indices = matrix->nonzero ? malloc(matrix->nonzero * sizeof(size_t)) : NULL;
  double *values = matrix->nonzero ? malloc(matrix->nonzero * sizeof(double)) : NULL;
  if (!offsets || !next || (matrix->nonzero && (!indices || !values))) {
    free(offsets); free(next); free(indices); free(values); return CRPX_ERROR_ALLOCATION;
  }
  for (size_t k = 0; k < matrix->nonzero; ++k) offsets[matrix->row_indices[k] + 1]++;
  for (size_t i = 0; i < matrix->rows; ++i) offsets[i + 1] += offsets[i];
  memcpy(next, offsets, matrix->rows * sizeof(size_t));
  for (size_t j = 0; j < matrix->columns; ++j)
    for (size_t k = matrix->offsets[j]; k < matrix->offsets[j + 1]; ++k) {
      size_t position = next[matrix->row_indices[k]]++;
      indices[position] = j;
      values[position] = matrix->values[k];
    }
  crpx_status status = crpx_sparse_matrix_new(matrix->columns, matrix->rows, matrix->nonzero,
                                             offsets, indices, values, out);
  free(offsets); free(next); free(indices); free(values);
  return status;
}

crpx_status crpx_sparse_matrix_multiply_vector(const struct crpx_sparse_matrix_struct *matrix, const double *input,
                                             size_t input_size, double *output, size_t output_size)
{
  if (!matrix || !input || !output) return CRPX_ERROR_INVALID_ARGUMENT;
  if (input_size != matrix->columns || output_size != matrix->rows) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < output_size; ++i) output[i] = 0;
  for (size_t j = 0; j < matrix->columns; ++j)
    for (size_t k = matrix->offsets[j]; k < matrix->offsets[j + 1]; ++k)
      output[matrix->row_indices[k]] += matrix->values[k] * input[j];
  return CRPX_OK;
}

void crpx_svd_free(crpx_svd_t decomposition)
{
  if (decomposition) { free(decomposition->singular); free(decomposition->left); free(decomposition->right); free(decomposition); }
}
size_t crpx_svd_rank(const struct crpx_svd_struct *decomposition) { return decomposition ? decomposition->rank : 0; }
size_t crpx_svd_rows(const struct crpx_svd_struct *decomposition) { return decomposition ? decomposition->rows : 0; }
size_t crpx_svd_columns(const struct crpx_svd_struct *decomposition) { return decomposition ? decomposition->columns : 0; }
const double *crpx_svd_singular_values(const struct crpx_svd_struct *decomposition)
{ return decomposition ? decomposition->singular : NULL; }
const double *crpx_svd_left_vectors(const struct crpx_svd_struct *decomposition)
{ return decomposition ? decomposition->left : NULL; }
const double *crpx_svd_right_vectors_transposed(const struct crpx_svd_struct *decomposition)
{ return decomposition ? decomposition->right : NULL; }

typedef struct { size_t index; double norm; } singular_index;
static int descending_norm(const void *left, const void *right)
{
  const singular_index *a = left, *b = right;
  if (a->norm != b->norm) return a->norm < b->norm ? 1 : -1;
  return (a->index > b->index) - (a->index < b->index);
}

crpx_status crpx_svd_decompose(const struct crpx_matrix_struct *matrix, size_t rank, double tolerance,
                              size_t sweeps, crpx_svd_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix || !isfinite(tolerance) || tolerance < 0 || tolerance >= 1) return CRPX_ERROR_INVALID_ARGUMENT;
  bool transposed = matrix->rows < matrix->columns;
  size_t m = transposed ? matrix->columns : matrix->rows, n = transposed ? matrix->rows : matrix->columns;
  if (!rank) rank = n;
  if (rank > n) return CRPX_ERROR_OUT_OF_RANGE;
  if (!product_fits(n, n, sizeof(double)) || n > SIZE_MAX / sizeof(singular_index)) return CRPX_ERROR_OVERFLOW;
  double scale = 0;
  for (size_t i = 0; i < matrix->rows * matrix->columns; ++i) {
    if (!isfinite(matrix->values[i])) return CRPX_ERROR_INVALID_ARGUMENT;
    scale = fmax(scale, fabs(matrix->values[i]));
  }
  if (scale == 0) scale = 1;
  double *work = malloc(m * n * sizeof(double)), *rotation = calloc(n * n, sizeof(double));
  singular_index *order = malloc(n * sizeof(*order));
  crpx_svd_t decomposition = calloc(1, sizeof(*decomposition));
  if (!work || !rotation || !order || !decomposition) {
    free(work); free(rotation); free(order); free(decomposition); return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      work[i * n + j] = (transposed ? matrix->values[j * matrix->columns + i] : matrix->values[i * n + j]) / scale;
  for (size_t j = 0; j < n; ++j) rotation[j * n + j] = 1;
  if (!tolerance) tolerance = 8 * DBL_EPSILON;
  if (!sweeps) sweeps = 100;
  bool converged = n == 1;
  for (size_t sweep = 0; sweep < sweeps; ++sweep) {
    bool rotated = false;
    for (size_t p = 0; p + 1 < n; ++p)
      for (size_t q = p + 1; q < n; ++q) {
        long double alpha = 0, beta = 0, gamma = 0;
        for (size_t i = 0; i < m; ++i) {
          long double a = work[i * n + p], b = work[i * n + q];
          alpha += a * a; beta += b * b; gamma += a * b;
        }
        if (gamma == 0 || fabsl(gamma) <= tolerance * sqrtl(alpha) * sqrtl(beta)) continue;
        long double tau = (beta - alpha) / (2 * gamma);
        long double tangent = copysignl(1, tau) / (fabsl(tau) + hypotl(1, tau));
        double cosine = (double)(1 / sqrtl(1 + tangent * tangent)), sine = (double)(tangent * cosine);
        for (size_t i = 0; i < m; ++i) {
          double a = work[i * n + p], b = work[i * n + q];
          work[i * n + p] = cosine * a - sine * b;
          work[i * n + q] = sine * a + cosine * b;
        }
        for (size_t i = 0; i < n; ++i) {
          double a = rotation[i * n + p], b = rotation[i * n + q];
          rotation[i * n + p] = cosine * a - sine * b;
          rotation[i * n + q] = sine * a + cosine * b;
        }
        rotated = true;
      }
    if (!rotated) { converged = true; break; }
  }
  if (!converged) {
    free(work); free(rotation); free(order); free(decomposition); return CRPX_ERROR_CONVERGENCE;
  }
  for (size_t j = 0; j < n; ++j) {
    long double sum = 0;
    for (size_t i = 0; i < m; ++i) sum += (long double)work[i * n + j] * work[i * n + j];
    order[j] = (singular_index){j, (double)sqrtl(sum)};
  }
  qsort(order, n, sizeof(*order), descending_norm);
  decomposition->rows = matrix->rows; decomposition->columns = matrix->columns; decomposition->rank = rank;
  decomposition->singular = malloc(rank * sizeof(double));
  decomposition->left = calloc(matrix->rows * rank, sizeof(double));
  decomposition->right = calloc(rank * matrix->columns, sizeof(double));
  if (!decomposition->singular || !decomposition->left || !decomposition->right) {
    free(work); free(rotation); free(order); crpx_svd_free(decomposition); return CRPX_ERROR_ALLOCATION;
  }
  crpx_status status = CRPX_OK;
  for (size_t k = 0; k < rank; ++k) {
    size_t j = order[k].index;
    double norm = order[k].norm;
    decomposition->singular[k] = norm * scale;
    if (!isfinite(decomposition->singular[k])) { status = CRPX_ERROR_OVERFLOW; break; }
    for (size_t i = 0; i < m; ++i) {
      double value = norm ? work[i * n + j] / norm : 0;
      if (transposed) decomposition->right[k * matrix->columns + i] = value;
      else decomposition->left[i * rank + k] = value;
    }
    for (size_t i = 0; i < n; ++i) {
      if (transposed) decomposition->left[i * rank + k] = rotation[i * n + j];
      else decomposition->right[k * matrix->columns + i] = rotation[i * n + j];
    }
  }
  free(work); free(rotation); free(order);
  if (status != CRPX_OK) { crpx_svd_free(decomposition); return status; }
  *out = decomposition;
  return CRPX_OK;
}

crpx_status crpx_svd_decompose_sparse(const struct crpx_sparse_matrix_struct *matrix, size_t rank,
                                     double tolerance, size_t sweeps, crpx_svd_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_matrix_t dense;
  crpx_status status = crpx_sparse_matrix_to_dense(matrix, &dense);
  if (status != CRPX_OK) return status;
  status = crpx_svd_decompose(dense, rank, tolerance, sweeps, out);
  crpx_matrix_free(dense);
  return status;
}

crpx_status crpx_svd_reconstruct(const struct crpx_svd_struct *decomposition, crpx_matrix_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!decomposition) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_matrix_new(decomposition->rows, decomposition->columns, out);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < decomposition->rows; ++i)
    for (size_t j = 0; j < decomposition->columns; ++j) {
      long double sum = 0;
      for (size_t k = 0; k < decomposition->rank; ++k)
        sum += (long double)decomposition->left[i * decomposition->rank + k] * decomposition->singular[k] *
               decomposition->right[k * decomposition->columns + j];
      (*out)->values[i * decomposition->columns + j] = (double)sum;
    }
  return CRPX_OK;
}
