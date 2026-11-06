/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "svd.h"
#include <math.h>
#include <stdio.h>
#include <stdint.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; } } while (0)

int main(void)
{
  for (size_t rows = 1; rows <= 12; ++rows)
    for (size_t columns = 1; columns <= 12; ++columns) {
      crpx_matrix_t matrix, reconstructed, transposed, roundtrip;
      CHECK(crpx_matrix_new(rows, columns, &matrix) == CRPX_OK);
      double *data = crpx_matrix_data(matrix);
      for (size_t i = 0; i < rows * columns; ++i) data[i] = sin((double)(i + 1) * 3.1) + i % 3;
      crpx_svd_t decomposition;
      CHECK(crpx_svd_decompose(matrix, 0, 0, 0, &decomposition) == CRPX_OK);
      CHECK(crpx_svd_reconstruct(decomposition, &reconstructed) == CRPX_OK);
      const double *copy = crpx_matrix_const_data(reconstructed);
      for (size_t i = 0; i < rows * columns; ++i) CHECK(fabs(copy[i] - data[i]) < 1e-12);
      const double *singular = crpx_svd_singular_values(decomposition);
      size_t rank = crpx_svd_rank(decomposition);
      for (size_t k = 1; k < rank; ++k) CHECK(singular[k - 1] >= singular[k]);
      const double *left = crpx_svd_left_vectors(decomposition);
      const double *right = crpx_svd_right_vectors_transposed(decomposition);
      for (size_t k = 0; k < rank; ++k)
        for (size_t l = 0; l < rank; ++l) {
          if (singular[k] < 1e-12 || singular[l] < 1e-12) continue;
          double u_dot = 0, v_dot = 0;
          for (size_t i = 0; i < rows; ++i) u_dot += left[i * rank + k] * left[i * rank + l];
          for (size_t i = 0; i < columns; ++i) v_dot += right[k * columns + i] * right[l * columns + i];
          CHECK(fabs(u_dot - (k == l)) < 1e-12);
          CHECK(fabs(v_dot - (k == l)) < 1e-12);
        }
      crpx_svd_free(decomposition);
      crpx_matrix_free(reconstructed);
      crpx_sparse_matrix_t sparse, sparse_transpose;
      CHECK(crpx_sparse_matrix_from_dense(matrix, 0, &sparse) == CRPX_OK);
      CHECK(crpx_sparse_matrix_transpose(sparse, &sparse_transpose) == CRPX_OK);
      CHECK(crpx_sparse_matrix_to_dense(sparse_transpose, &transposed) == CRPX_OK);
      CHECK(crpx_matrix_transpose(transposed, &roundtrip) == CRPX_OK);
      for (size_t i = 0; i < rows * columns; ++i) CHECK(crpx_matrix_const_data(roundtrip)[i] == data[i]);
      double input[12], dense_output[12], sparse_output[12];
      for (size_t i = 0; i < columns; ++i) input[i] = (double)i - 4;
      CHECK(crpx_matrix_multiply_vector(matrix, input, columns, dense_output, rows) == CRPX_OK);
      CHECK(crpx_sparse_matrix_multiply_vector(sparse, input, columns, sparse_output, rows) == CRPX_OK);
      for (size_t i = 0; i < rows; ++i) CHECK(fabs(dense_output[i] - sparse_output[i]) < 1e-12);
      crpx_sparse_matrix_free(sparse); crpx_sparse_matrix_free(sparse_transpose);
      crpx_matrix_free(transposed); crpx_matrix_free(roundtrip); crpx_matrix_free(matrix);
    }
  double diagonal[] = {4, 0, 0, 0, 3, 0};
  crpx_matrix_t matrix, reconstructed;
  crpx_svd_t decomposition;
  CHECK(crpx_matrix_from_array(diagonal, 2, 3, &matrix) == CRPX_OK);
  CHECK(crpx_svd_decompose(matrix, 1, 0, 0, &decomposition) == CRPX_OK);
  CHECK(fabs(crpx_svd_singular_values(decomposition)[0] - 4) < 1e-15);
  CHECK(crpx_svd_reconstruct(decomposition, &reconstructed) == CRPX_OK);
  CHECK(crpx_matrix_const_data(reconstructed)[0] == 4);
  CHECK(crpx_matrix_const_data(reconstructed)[4] == 0);
  crpx_matrix_free(reconstructed); crpx_svd_free(decomposition); crpx_matrix_free(matrix);
  CHECK(crpx_matrix_new(4, 4, &matrix) == CRPX_OK);
  for (size_t i = 0; i < 16; ++i) crpx_matrix_data(matrix)[i] = (double)(i + 1);
  CHECK(crpx_svd_decompose(matrix, 0, 0, 1, &decomposition) == CRPX_ERROR_CONVERGENCE && !decomposition);
  CHECK(crpx_svd_decompose(matrix, 0, 0, 0, &decomposition) == CRPX_OK);
  CHECK(crpx_svd_reconstruct(decomposition, &reconstructed) == CRPX_OK);
  for (size_t i = 0; i < 16; ++i) CHECK(fabs(crpx_matrix_const_data(reconstructed)[i] - (double)(i + 1)) < 1e-12);
  crpx_svd_free(decomposition); crpx_matrix_free(reconstructed);
  crpx_matrix_data(matrix)[0] = NAN;
  CHECK(crpx_svd_decompose(matrix, 0, 0, 0, &decomposition) == CRPX_ERROR_INVALID_ARGUMENT && !decomposition);
  crpx_matrix_free(matrix);
  CHECK(crpx_matrix_new(SIZE_MAX, 2, &matrix) == CRPX_ERROR_OVERFLOW && !matrix);
  size_t offsets[] = {0, 2}, rows[] = {1, 1}; double values[] = {1, 2};
  crpx_sparse_matrix_t sparse;
  CHECK(crpx_sparse_matrix_new(2, 1, 2, offsets, rows, values, &sparse) == CRPX_ERROR_INVALID_ARGUMENT && !sparse);
  return 0;
}
