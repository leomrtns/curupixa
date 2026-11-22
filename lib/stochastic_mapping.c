/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
#include "stochastic_mapping.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
struct crpx_mapping_workspace {
  double *matrices, *backward, *weights, *state_weights;
  size_t matrix_capacity, backward_capacity, weight_capacity, state_capacity;
};
static crpx_status reserve (double **array, size_t *capacity, size_t count);
static int compare_times (const void *a, const void *b);
crpx_status
crpx_mapping_workspace_new (crpx_mapping_workspace **out)
{
  /* Choose the virtual-event count conditional on endpoints, then sample its embedded bridge. */
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc (1, sizeof (**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void
crpx_mapping_workspace_free (crpx_mapping_workspace *workspace)
{
  if (workspace) {
    free (workspace->matrices);
    free (workspace->backward);
    free (workspace->weights);
    free (workspace->state_weights);
    free (workspace);
  }
}
void
crpx_branch_history_free (crpx_branch_history *history)
{
  if (history) {
    free (history->times);
    free (history->states);
    free (history);
  }
}
crpx_status
crpx_substitution_bridge_sample (const crpx_substitution_model *model, double length, size_t first, size_t last,
                                 double tolerance, size_t max_virtual_events, crpx_uniform_function uniform,
                                 void *state, crpx_mapping_workspace *workspace, crpx_branch_history **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t states = crpx_substitution_model_states (model);
  if (!model || !workspace || !uniform || !isfinite (length) || length < 0 || first >= states || last >= states ||
      !isfinite (tolerance) || tolerance <= 0 || tolerance >= 1)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (max_virtual_events == SIZE_MAX || states > SIZE_MAX / states / 2 / sizeof (double) ||
      max_virtual_events + 1 > SIZE_MAX / states / sizeof (double))
    return CRPX_ERROR_OVERFLOW;
  crpx_status status = reserve (&workspace->matrices, &workspace->matrix_capacity, 2 * states * states);
  if (status == CRPX_OK)
    status = reserve (&workspace->backward, &workspace->backward_capacity, (max_virtual_events + 1) * states);
  if (status == CRPX_OK) status = reserve (&workspace->weights, &workspace->weight_capacity, max_virtual_events + 1);
  if (status == CRPX_OK) status = reserve (&workspace->state_weights, &workspace->state_capacity, states);
  if (status != CRPX_OK) return status;
  double *generator = workspace->matrices, *transition = generator + states * states;
  status = crpx_substitution_model_generator (model, generator, states * states);
  if (status == CRPX_OK)
    status = crpx_substitution_model_transition (model, length, CRPX_BRANCH_FIXED, transition, states * states);
  if (status != CRPX_OK) return status;
  double probability = transition[first * states + last], rate = 0;
  if (probability <= 0) return CRPX_ERROR_NOT_FOUND;
  for (size_t i = 0; i < states; i++)
    if (-generator[i * states + i] > rate) rate = -generator[i * states + i];
  double mean = rate * length;
  if (!isfinite (mean)) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < states; i++)
    for (size_t j = 0; j < states; j++)
      transition[i * states + j] = rate ? generator[i * states + j] / rate + (i == j) : (i == j);
  for (size_t i = 0; i < states; i++)
    workspace->backward[i] = i == last;
  long double cumulative = 0;
  size_t terms = 0;
  double bound = INFINITY;
  for (size_t events = 0; events <= max_virtual_events; events++) {
    if (events)
      for (size_t i = 0; i < states; i++) {
        long double value = 0;
        for (size_t j = 0; j < states; j++)
          value += (long double) transition[i * states + j] * workspace->backward[(events - 1) * states + j];
        workspace->backward[events * states + i] = (double) value;
      }
    double log_poisson = mean ? -mean + events * log (mean) - lgamma ((double) events + 1) : (events ? -INFINITY : 0);
    workspace->weights[events] = log_poisson + log (workspace->backward[events * states + first]);
    cumulative += expl (log_poisson);
    bound = (double) ((fmaxl (0, 1 - cumulative) + 64 * DBL_EPSILON * (events + 1)) / probability);
    terms = events + 1;
    if ((double) events >= mean && bound <= tolerance) break;
  }
  if (bound > tolerance) return CRPX_ERROR_CONVERGENCE;
  size_t events;
  status = crpx_sample_log_weights (workspace->weights, terms, uniform, state, &events);
  if (status != CRPX_OK) return status;
  crpx_branch_history *history = calloc (1, sizeof (*history));
  if (!history) return CRPX_ERROR_ALLOCATION;
  history->times = events ? malloc (events * sizeof (double)) : NULL;
  history->states = malloc ((events + 1) * sizeof (size_t));
  if ((events && !history->times) || !history->states) {
    status = CRPX_ERROR_ALLOCATION;
    goto done;
  }
  for (size_t i = 0; i < events; i++) {
    double draw;
    status = crpx_sample_uniform (uniform, state, &draw);
    if (status != CRPX_OK) goto done;
    /* Zero has measure zero under a continuous uniform; map it to an interior double. */
    history->times[i] = length * (draw ? draw : 0x1p-54);
  }
  if (events) qsort (history->times, events, sizeof (double), compare_times);
  size_t current = first;
  history->states[0] = current;
  history->omitted_probability_bound = bound;
  for (size_t i = 0; i < events; i++) {
    for (size_t next = 0; next < states; next++)
      workspace->state_weights[next] =
        log (transition[current * states + next]) + log (workspace->backward[(events - i - 1) * states + next]);
    size_t next;
    status = crpx_sample_log_weights (workspace->state_weights, states, uniform, state, &next);
    if (status != CRPX_OK) goto done;
    if (next != current) {
      history->times[history->count] = history->times[i];
      history->states[++history->count] = next;
    }
    current = next;
  }
  *out = history;
  return CRPX_OK;
done:
  crpx_branch_history_free (history);
  return status;
}
static crpx_status
reserve (double **array, size_t *capacity, size_t count)
{
  if (count <= *capacity) return CRPX_OK;
  if (count > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  double *replacement = realloc (*array, count * sizeof (double));
  if (!replacement) return CRPX_ERROR_ALLOCATION;
  *array = replacement;
  *capacity = count;
  return CRPX_OK;
}
static int
compare_times (const void *a, const void *b)
{
  double first = *(const double *) a, second = *(const double *) b;
  return (first > second) - (first < second);
}
