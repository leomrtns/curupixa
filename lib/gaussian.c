/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "gaussian.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
struct crpx_regression_workspace {
  double *matrix, *response, *inverse;
  size_t cells, rows, square;
};
static crpx_status reserve (double **array, size_t *capacity, size_t required);
void
crpx_covariance_free (crpx_covariance *covariance)
{
  if (covariance) {
    free (covariance->lower);
    free (covariance);
  }
}
crpx_status
crpx_covariance_new (const double *matrix, size_t size, crpx_covariance **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!matrix || !size) return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / size / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  crpx_covariance *factor = calloc (1, sizeof (*factor));
  if (!factor) return CRPX_ERROR_ALLOCATION;
  factor->size = size;
  factor->lower = calloc (size * size, sizeof (double));
  if (!factor->lower) {
    crpx_covariance_free (factor);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i = 0; i < size; i++)
    for (size_t j = 0; j <= i; j++) {
      if (!isfinite (matrix[i * size + j]) || !isfinite (matrix[j * size + i]) ||
          fabs (matrix[i * size + j] - matrix[j * size + i]) >
            1e-12 * fmax (1, fmax (fabs (matrix[i * size + j]), fabs (matrix[j * size + i])))) {
        crpx_covariance_free (factor);
        return CRPX_ERROR_INVALID_ARGUMENT;
      }
      long double sum = matrix[i * size + j];
      for (size_t k = 0; k < j; k++)
        sum -= (long double) factor->lower[i * size + k] * factor->lower[j * size + k];
      if (i == j) {
        if (sum <= 0 || !isfinite (sum)) {
          crpx_covariance_free (factor);
          return CRPX_ERROR_INVALID_ARGUMENT;
        }
        factor->lower[i * size + j] = (double) sqrtl (sum);
        factor->log_determinant += log ((double) sum);
      } else
        factor->lower[i * size + j] = (double) (sum / factor->lower[j * size + j]);
    }
  *out = factor;
  return CRPX_OK;
}
crpx_status
crpx_regression_workspace_new (crpx_regression_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void
crpx_regression_workspace_free (crpx_regression_workspace *workspace)
{
  if (workspace) {
    free (workspace->matrix);
    free (workspace->response);
    free (workspace->inverse);
    free (workspace);
  }
}
void
crpx_regression_result_free (crpx_regression_result *result)
{
  if (result) {
    free (result->coefficients);
    free (result->coefficient_covariance);
    free (result);
  }
}
crpx_status
crpx_gaussian_regression (const crpx_covariance *factor, const double *design, size_t predictors,
                          const double *response, bool reml, crpx_regression_workspace *workspace,
                          crpx_regression_result **out)
{ /* Whiten by the covariance factor, then solve and estimate uncertainty by QR. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!factor || !design || !response || !workspace || !predictors || predictors >= factor->size)
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t rows = factor->size, columns = predictors;
  if (rows > SIZE_MAX / columns / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = reserve (&workspace->matrix, &workspace->cells, rows * columns);
  if (status == CRPX_OK) status = reserve (&workspace->response, &workspace->rows, rows);
  if (status == CRPX_OK) status = reserve (&workspace->inverse, &workspace->square, columns * columns);
  if (status != CRPX_OK) return status;
  double *matrix = workspace->matrix, *values = workspace->response, *inverse = workspace->inverse;
  for (size_t i = 0; i < rows; i++) {
    if (!isfinite (response[i])) return CRPX_ERROR_INVALID_ARGUMENT;
    long double value = response[i];
    for (size_t k = 0; k < i; k++)
      value -= (long double) factor->lower[i * rows + k] * values[k];
    values[i] = (double) (value / factor->lower[i * rows + i]);
    if (!isfinite (values[i])) return CRPX_ERROR_OVERFLOW;
    for (size_t j = 0; j < columns; j++) {
      if (!isfinite (design[i * columns + j])) return CRPX_ERROR_INVALID_ARGUMENT;
      long double value = design[i * columns + j];
      for (size_t k = 0; k < i; k++)
        value -= (long double) factor->lower[i * rows + k] * matrix[k * columns + j];
      matrix[i * columns + j] = (double) (value / factor->lower[i * rows + i]);
    }
  }
  double design_scale = 0;
  for (size_t i = 0; i < rows * columns; i++) {
    if (!isfinite (matrix[i])) return CRPX_ERROR_OVERFLOW;
    if (fabs (matrix[i]) > design_scale) design_scale = fabs (matrix[i]);
  }
  double log_design = 0;
  for (size_t column = 0; column < columns; column++) {
    double norm = 0;
    for (size_t i = column; i < rows; i++)
      norm = hypot (norm, matrix[i * columns + column]);
    if (!isfinite (norm) || norm <= 64 * DBL_EPSILON * design_scale * sqrt ((double) rows))
      return CRPX_ERROR_INVALID_ARGUMENT;
    double old = matrix[column * columns + column], diagonal = -copysign (norm, old);
    double first = old - diagonal, tau = (diagonal - old) / diagonal;
    if (!isfinite (first)) return CRPX_ERROR_OVERFLOW;
    matrix[column * columns + column] = diagonal;
    for (size_t i = column + 1; i < rows; i++)
      matrix[i * columns + column] /= first;
    for (size_t j = column + 1; j < columns; j++) {
      long double product = matrix[column * columns + j];
      for (size_t i = column + 1; i < rows; i++)
        product += (long double) matrix[i * columns + column] * matrix[i * columns + j];
      double change = (double) (tau * product);
      matrix[column * columns + j] -= change;
      for (size_t i = column + 1; i < rows; i++)
        matrix[i * columns + j] -= matrix[i * columns + column] * change;
    }
    long double product = values[column];
    for (size_t i = column + 1; i < rows; i++)
      product += (long double) matrix[i * columns + column] * values[i];
    double change = (double) (tau * product);
    values[column] -= change;
    for (size_t i = column + 1; i < rows; i++)
      values[i] -= matrix[i * columns + column] * change;
    log_design += 2 * log (fabs (diagonal));
  }
  crpx_regression_result *result = calloc (1, sizeof (*result));
  if (!result) return CRPX_ERROR_ALLOCATION;
  result->observations = rows;
  result->predictors = columns;
  result->coefficients = malloc (columns * sizeof (double));
  result->coefficient_covariance = calloc (columns * columns, sizeof (double));
  if (!result->coefficients || !result->coefficient_covariance) {
    crpx_regression_result_free (result);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t position = columns; position > 0;) {
    size_t i = --position;
    long double value = values[i];
    for (size_t j = i + 1; j < columns; j++)
      value -= (long double) matrix[i * columns + j] * result->coefficients[j];
    result->coefficients[i] = (double) (value / matrix[i * columns + i]);
  }
  long double residual = 0;
  for (size_t i = columns; i < rows; i++)
    residual += (long double) values[i] * values[i];
  result->residual_sum_squares = (double) residual;
  double degrees = (double) (reml ? rows - columns : rows);
  result->variance = (double) (residual / degrees);
  result->log_likelihood = residual == 0 ? INFINITY
                                         : -0.5 * (degrees * (log (2 * acos (-1.0) * result->variance) + 1) +
                                                   factor->log_determinant + (reml ? log_design : 0));
  memset (inverse, 0, columns * columns * sizeof (double));
  for (size_t column = 0; column < columns; column++)
    for (size_t position = columns; position > 0;) {
      size_t row = --position;
      long double value = row == column ? 1 : 0;
      for (size_t j = row + 1; j < columns; j++)
        value -= (long double) matrix[row * columns + j] * inverse[j * columns + column];
      inverse[row * columns + column] = (double) (value / matrix[row * columns + row]);
    }
  for (size_t i = 0; i < columns; i++)
    for (size_t j = 0; j < columns; j++) {
      long double value = 0;
      for (size_t k = 0; k < columns; k++)
        value += (long double) inverse[i * columns + k] * inverse[j * columns + k];
      result->coefficient_covariance[i * columns + j] = (double) (value * result->variance);
    }
  if (!isfinite (result->residual_sum_squares) || !isfinite (result->variance) ||
      (residual > 0 && !isfinite (result->log_likelihood))) {
    crpx_regression_result_free (result);
    return CRPX_ERROR_OVERFLOW;
  }
  for (size_t i = 0; i < columns; i++)
    if (!isfinite (result->coefficients[i])) {
      crpx_regression_result_free (result);
      return CRPX_ERROR_OVERFLOW;
    }
  *out = result;
  return CRPX_OK;
}

crpx_status
crpx_covariance_solve (const crpx_covariance *factor, const double *input, double *output, size_t count)
{
  if (!factor || !input || !output || count != factor->size) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < count; i++) {
    long double value = input[i];
    if (!isfinite (value)) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t j = 0; j < i; j++)
      value -= (long double) factor->lower[i * count + j] * output[j];
    output[i] = (double) (value / factor->lower[i * count + i]);
  }
  for (size_t position = count; position > 0;) {
    size_t i = --position;
    long double value = output[i];
    for (size_t j = i + 1; j < count; j++)
      value -= (long double) factor->lower[j * count + i] * output[j];
    output[i] = (double) (value / factor->lower[i * count + i]);
    if (!isfinite (output[i])) return CRPX_ERROR_OVERFLOW;
  }
  return CRPX_OK;
}
void
crpx_logistic_result_free (crpx_logistic_result *result)
{
  if (result) {
    free (result->coefficients);
    free (result->latent);
    free (result->probabilities);
    free (result);
  }
}
crpx_status
crpx_logistic_mixed_fit (const crpx_covariance *random_covariance, const double *design, size_t predictors,
                         const uint8_t *response, double tolerance, size_t max_iterations,
                         crpx_regression_workspace *workspace, crpx_logistic_result **out)
{ /* PQL alternates working Gaussian fits with latent Gaussian predictions. The covariance is supplied, not fitted. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!random_covariance || !design || !response || !workspace || !predictors ||
      predictors >= random_covariance->size || !isfinite (tolerance) || tolerance <= 0 || !max_iterations)
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t rows = random_covariance->size;
  if (rows > SIZE_MAX / rows / sizeof (double) || rows > SIZE_MAX / 6 / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double *matrix = malloc (rows * rows * sizeof (double));
  double *buffer = calloc (6 * rows, sizeof (double));
  crpx_logistic_result *result = calloc (1, sizeof (*result));
  crpx_status status = CRPX_OK;
  if (!matrix || !buffer || !result) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  result->observations = rows;
  result->predictors = predictors;
  result->coefficients = calloc (predictors, sizeof (double));
  result->latent = calloc (rows, sizeof (double));
  result->probabilities = calloc (rows, sizeof (double));
  if (!result->coefficients || !result->latent || !result->probabilities) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  double *eta = buffer, *working = eta + rows, *weights = working + rows, *residual = weights + rows;
  double *solution = residual + rows, *old = solution + rows;
  for (size_t i = 0; i < rows; i++)
    if (response[i] > 1) {
      status = CRPX_ERROR_INVALID_ARGUMENT;
      goto done;
    }
  bool converged = false;
  for (size_t iteration = 0; iteration < max_iterations; iteration++) {
    for (size_t i = 0; i < rows; i++) {
      double probability = eta[i] >= 0 ? 1 / (1 + exp (-eta[i])) : exp (eta[i]) / (1 + exp (eta[i]));
      probability = fmin (1 - 1e-12, fmax (1e-12, probability));
      weights[i] = probability * (1 - probability);
      working[i] = eta[i] + (response[i] - probability) / weights[i];
      old[i] = eta[i];
      for (size_t j = 0; j <= i; j++) {
        long double value = 0;
        for (size_t k = 0; k <= j; k++)
          value += (long double) random_covariance->lower[i * rows + k] * random_covariance->lower[j * rows + k];
        matrix[i * rows + j] = matrix[j * rows + i] = (double) value;
      }
      matrix[i * rows + i] += 1 / weights[i];
    }
    crpx_covariance *working_covariance = NULL;
    crpx_regression_result *fit = NULL;
    status = crpx_covariance_new (matrix, rows, &working_covariance);
    if (status == CRPX_OK)
      status = crpx_gaussian_regression (working_covariance, design, predictors, working, false, workspace, &fit);
    if (status != CRPX_OK) {
      crpx_covariance_free (working_covariance);
      crpx_regression_result_free (fit);
      goto done;
    }
    memcpy (result->coefficients, fit->coefficients, predictors * sizeof (double));
    for (size_t i = 0; i < rows; i++) {
      residual[i] = working[i];
      for (size_t j = 0; j < predictors; j++)
        residual[i] -= design[i * predictors + j] * result->coefficients[j];
    }
    status = crpx_covariance_solve (working_covariance, residual, solution, rows);
    crpx_covariance_free (working_covariance);
    crpx_regression_result_free (fit);
    if (status != CRPX_OK) goto done;
    double maximum_change = 0;
    for (size_t i = 0; i < rows; i++) {
      long double latent = 0;
      for (size_t j = 0; j < rows; j++)
        latent += (long double) (matrix[i * rows + j] - (i == j ? 1 / weights[i] : 0)) * solution[j];
      result->latent[i] = (double) latent;
      eta[i] = result->latent[i];
      for (size_t j = 0; j < predictors; j++)
        eta[i] += design[i * predictors + j] * result->coefficients[j];
      if (!isfinite (eta[i])) {
        status = CRPX_ERROR_OVERFLOW;
        goto done;
      }
      double change = fabs (eta[i] - old[i]) / fmax (1, fabs (old[i]));
      if (change > maximum_change) maximum_change = change;
    }
    result->iterations = iteration + 1;
    if (maximum_change < tolerance) {
      converged = true;
      break;
    }
  }
  if (!converged) {
    status = CRPX_ERROR_CONVERGENCE;
    goto done;
  }
  status = crpx_covariance_solve (random_covariance, result->latent, solution, rows);
  if (status != CRPX_OK) goto done;
  long double log_density = 0;
  for (size_t i = 0; i < rows; i++) {
    result->probabilities[i] = eta[i] >= 0 ? 1 / (1 + exp (-eta[i])) : exp (eta[i]) / (1 + exp (eta[i]));
    log_density +=
      response[i] * eta[i] - fmax (eta[i], 0) - log1p (exp (-fabs (eta[i]))) - 0.5L * result->latent[i] * solution[i];
  }
  result->log_penalised_likelihood = (double) log_density;
  *out = result;
  result = NULL;
done:
  free (matrix);
  free (buffer);
  crpx_logistic_result_free (result);
  return status;
}

void
crpx_gaussian_prediction_free (crpx_gaussian_prediction *prediction)
{
  if (prediction) {
    free (prediction->gain);
    free (prediction->variances);
    free (prediction);
  }
}
crpx_status
crpx_gaussian_prediction_new (const crpx_covariance *observed, const double *cross_covariance,
                              const double *target_variances, size_t targets, crpx_gaussian_prediction **out)
{ /* Condition a Gaussian vector once; later response vectors reuse the same gains and variances. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!observed || !cross_covariance || !target_variances || !targets) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t count = observed->size;
  if (targets > SIZE_MAX / count / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  crpx_gaussian_prediction *prediction = calloc (1, sizeof (*prediction));
  if (!prediction) return CRPX_ERROR_ALLOCATION;
  prediction->observations = count;
  prediction->targets = targets;
  prediction->gain = malloc (targets * count * sizeof (double));
  prediction->variances = malloc (targets * sizeof (double));
  if (!prediction->gain || !prediction->variances) {
    crpx_gaussian_prediction_free (prediction);
    return CRPX_ERROR_ALLOCATION;
  }
  for (size_t target = 0; target < targets; target++) {
    if (!isfinite (target_variances[target]) || target_variances[target] < 0) {
      crpx_gaussian_prediction_free (prediction);
      return CRPX_ERROR_INVALID_ARGUMENT;
    }
    crpx_status status =
      crpx_covariance_solve (observed, cross_covariance + target * count, prediction->gain + target * count, count);
    if (status != CRPX_OK) {
      crpx_gaussian_prediction_free (prediction);
      return status;
    }
    long double variance = target_variances[target];
    for (size_t i = 0; i < count; i++)
      variance -= (long double) prediction->gain[target * count + i] * cross_covariance[target * count + i];
    if (!isfinite (variance) || variance < -1e-10 * fmax (1, target_variances[target])) {
      crpx_gaussian_prediction_free (prediction);
      return CRPX_ERROR_INVALID_ARGUMENT;
    }
    prediction->variances[target] = (double) fmaxl (0, variance);
  }
  *out = prediction;
  return CRPX_OK;
}
crpx_status
crpx_gaussian_predict (const crpx_gaussian_prediction *prediction, const double *response, const double *observed_means,
                       const double *target_means, double *means, double *variances)
{
  if (!prediction || !response || !observed_means || !target_means || !means || !variances)
    return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < prediction->observations; i++)
    if (!isfinite (response[i]) || !isfinite (observed_means[i])) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t target = 0; target < prediction->targets; target++)
    if (!isfinite (target_means[target])) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t target = 0; target < prediction->targets; target++) {
    long double value = target_means[target];
    for (size_t i = 0; i < prediction->observations; i++)
      value +=
        (long double) prediction->gain[target * prediction->observations + i] * (response[i] - observed_means[i]);
    if (!isfinite ((double) value)) return CRPX_ERROR_OVERFLOW;
    means[target] = (double) value;
    variances[target] = prediction->variances[target];
  }
  return CRPX_OK;
}

static crpx_status
reserve (double **array, size_t *capacity, size_t required)
{
  if (required <= *capacity) return CRPX_OK;
  if (required > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double *replacement = realloc (*array, required * sizeof (double));
  if (!replacement) return CRPX_ERROR_ALLOCATION;
  *array = replacement;
  *capacity = required;
  return CRPX_OK;
}
