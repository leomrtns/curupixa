/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file parsimony_internal.h
    \brief Shared state-cost transformations for scoring and ancestral reconstruction.
    Sankoff (1975), doi:10.1137/0128004. Internal arrays must contain states entries. */
#ifndef CURUPIXA_PARSIMONY_INTERNAL_H
#define CURUPIXA_PARSIMONY_INTERNAL_H
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
/* mode: 0 unordered, 1 ordered absolute difference, 2 custom directed matrix. */
static inline double
crpx_private_parsimony_add (double first, double second)
{
  long double sum = (long double) first + second;
  return sum > DBL_MAX ? INFINITY : (double) sum;
}
static inline void
crpx_private_parsimony_transform (const double *input, double *output, size_t states, unsigned mode,
                                  const double *costs, bool transpose)
{
  if (mode == 0) {
    double minimum = INFINITY;
    for (size_t state = 0; state < states; state++)
      if (input[state] < minimum) minimum = input[state];
    for (size_t state = 0; state < states; state++)
      output[state] = fmin (input[state], crpx_private_parsimony_add (minimum, 1));
  } else if (mode == 1) {
    output[0] = input[0];
    for (size_t state = 1; state < states; state++)
      output[state] = fmin (input[state], crpx_private_parsimony_add (output[state - 1], 1));
    for (size_t state = states - 1; state > 0; state--)
      output[state - 1] = fmin (output[state - 1], crpx_private_parsimony_add (output[state], 1));
  } else
    for (size_t state = 0; state < states; state++) {
      output[state] = INFINITY;
      for (size_t next = 0; next < states; next++) {
        double cost = costs[transpose ? next * states + state : state * states + next];
        double value = crpx_private_parsimony_add (input[next], cost);
        if (value < output[state]) output[state] = value;
      }
    }
}
/* Integer specialisation keeps large weighted mask scores exact without converting to double. */
static inline void
crpx_private_parsimony_unit_accumulate (uint64_t *parent, const uint64_t *child, size_t states, uint64_t infinity)
{
  uint64_t minimum = infinity;
  for (size_t state = 0; state < states; state++)
    if (child[state] < minimum) minimum = child[state];
  for (size_t state = 0; state < states; state++) {
    uint64_t value = child[state] < minimum + 1 ? child[state] : minimum + 1;
    parent[state] = parent[state] >= infinity - value ? infinity : parent[state] + value;
  }
}
#endif
