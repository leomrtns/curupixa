/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/* Compatibility interface over the checked, seekable numerical sequence module. */
#include "quasi_random.h"
#include <math.h>
crpx_quasi_random_t
new_crpx_quasi_random (crpx_global_t context, size_t dimensions)
{
  if (!context || !dimensions || dimensions > 2048) return NULL;
  crpx_quasi_random_t q = calloc (1, sizeof (*q));
  if (!q) return NULL;
  q->size = dimensions;
  q->rem = dimensions;
  q->cglob = crpx_global_retain (context);
  q->r = calloc (dimensions, sizeof (double));
  q->ko = malloc (dimensions * sizeof (double));
  if (!q->r || !q->ko || crpx_quasi_sequence_new (CRPX_QUASI_HALTON, dimensions, &q->halton) != CRPX_OK ||
      crpx_quasi_sequence_new (CRPX_QUASI_REVERSE_HALTON, dimensions, &q->reverse) != CRPX_OK) {
    del_crpx_quasi_random (q);
    return NULL;
  }
  crpx_quasi_random_reset (q);
  return q;
}
void
del_crpx_quasi_random (crpx_quasi_random_t q)
{
  if (!q) return;
  crpx_quasi_sequence_free (q->halton);
  crpx_quasi_sequence_free (q->reverse);
  free (q->r);
  free (q->ko);
  crpx_global_finalise (q->cglob);
  free (q);
}
void
crpx_quasi_random_reset (crpx_quasi_random_t q)
{
  if (!q) return;
  q->iteration = 0;
  for (size_t i = 0; i < q->size; i++) {
    q->r[i] = 0;
    q->ko[i] = crpx_random_double_positive (q->cglob);
  }
}
void
crpx_quasi_random_next_korobov (crpx_quasi_random_t q)
{
  if (!q) return;
  q->iteration++;
  for (size_t i = 0; i < q->size; i++) {
    q->r[i] += q->ko[i];
    if (q->r[i] >= 1) q->r[i] -= 1;
  }
}
void
crpx_quasi_random_next_halton (crpx_quasi_random_t q)
{
  if (!q) return;
  crpx_status status = crpx_quasi_sequence_seek (q->reverse, ++q->iteration);
  if (status == CRPX_OK) status = crpx_quasi_sequence_next (q->reverse, q->r, q->size);
  if (status != CRPX_OK)
    for (size_t i = 0; i < q->size; i++)
      q->r[i] = NAN;
}
void
crpx_quasi_random_next_halton_original (crpx_quasi_random_t q)
{
  if (!q) return;
  crpx_status status = crpx_quasi_sequence_seek (q->halton, ++q->iteration);
  if (status == CRPX_OK) status = crpx_quasi_sequence_next (q->halton, q->r, q->size);
  if (status != CRPX_OK)
    for (size_t i = 0; i < q->size; i++)
      q->r[i] = NAN;
  /* Preserve the original interface's reversed coordinate order. */
  for (size_t i = 0; i < q->size / 2; i++) {
    double x = q->r[i];
    q->r[i] = q->r[q->size - i - 1];
    q->r[q->size - i - 1] = x;
  }
}
