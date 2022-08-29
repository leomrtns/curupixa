/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file quasi_random.c 
 *  \brief Quasi-random number generator. Based on algorithms from the Gnu Scientific Library (GPL-3) and R's qrng*/

#include "quasi_random.h"
#include "quasi_random_constants.h"

/* Copyright (C) 2007 O. Teytaud olivier.teytaud@inria.fr for the GSL under the GPL-3 licence 
 * Implementation for Halton generator, modified according to the following improvement (reverse scrambling): [ B. Vandewoestyne and
 * R. Cools, Good permutations for deterministic scrambled Halton sequences in terms of L2-discrepancy, Computational and Applied
 * Mathematics 189, 2006]  */
static double reverse_vdcorput (int x, int b);
static double original_vdcorput (int x, int b); // original, without the reverse scrambling

crpx_quasi_random_t
new_crpx_quasi_random (crpx_global_t cglob, size_t size)
{
  crpx_quasi_random_t q = (crpx_quasi_random_t) crpx_malloc (cglob, sizeof (crpx_quasi_random_struct));
  q->r = q->ko = NULL;
  q->r = (double*) crpx_malloc (cglob, size  * sizeof (double));
  if (!q->r) {
    free (q);
    return NULL;
  }
  q->ko = (double *) crpx_malloc (cglob, size  * sizeof (double));
  if (!q->ko) {
    free (q);
    return NULL;
  }
  q->size = size;
  q->iteration = 0;
  q->leap = size / HALTON_MAX_DIMENSION;
  q->rem  = size % HALTON_MAX_DIMENSION;

  if (size > HALTON_MAX_DIMENSION) crpx_logger_warning (cglob, "Halton quasi-random generator is not efficient for space dimensions > %d", HALTON_MAX_DIMENSION);
  crpx_link_add_global_pointer (cglob, q->cglob); // thread-safe increase of ref_counter
  crpx_quasi_random_reset (q);
  return q;
}

void
del_crpx_quasi_random (crpx_quasi_random_t q)
{
  if (!q) return;
  crpx_global_finalise (q->cglob); // it should just decrease ref_counter unless user called it in wrong order 
  if (q->r)  free (q->r);
  if (q->ko) free (q->ko);
  free (q);
}

void
crpx_quasi_random_reset (crpx_quasi_random_t q)
{
  if (!q) return;
  for (size_t i = 0; i < q->size; i++) {
    q->r[i] = 0.;
    q->ko[i] = crpx_random_double_positive (q->cglob); // original <github.com/cran/qrng> uses user-defined generator smaller than total iterations
  }
  q->iteration = 0;
  return;
}

void
crpx_quasi_random_next_korobov (crpx_quasi_random_t q)
{ // https://github.com/cran/qrng/blob/master/src/korobov.c
  if (!q) return;
  q->iteration++;
  for (size_t i = 0; i < q->size; i++) {
    q->r[i] += q->ko[i];
    if (q->r[i] >= 1.) q->r[i] -= 1.;
  }
  return;
}

/*! \brief Returns n_qrand quasi-random numbers in [0,1), skipping trivial [0,0,..,0]; iteration should be initialised to zero; originally works only if n_qrand < MAX_DIMENSION = 1229 
 *  here we add more dimensions by leaping */
void
crpx_quasi_random_next_halton (crpx_quasi_random_t q)
{
  size_t i, j, k=0, iter_skip;
  q->iteration++;
  iter_skip = q->iteration;

  for (k = 0, i = 0; i < q->leap; i++) { // only enter this loop if q->size > HALTON_MAX_DIMENSION
    if (i < HALTON_MAX_DIMENSION + 1) iter_skip = (q->iteration) * halton_prime_numbers[i] + i; // "x prime" means "leap" and "+i" means "warm up"
    else iter_skip = (q->iteration) * i; // we're too far the original assumption of 1229 dimensions
    for (j = 1; j <= HALTON_MAX_DIMENSION; j++) {
      q->r[k++] = reverse_vdcorput (iter_skip, halton_prime_numbers[j]); // original GSL has prime_number[0...h-1] instead of [1...h] but we use prime_number[0] to store "1"
    }
  }

  for (j = 1; j <= q->rem; j++) q->r[k++] = reverse_vdcorput (iter_skip, halton_prime_numbers[j]); // original GSL has prime_number[j] instead of j+1
  assert ((k == (q->size-1)) && "not all quasi-random elements were generated in crpx_quasi_random_next_halton");
  return;
}

void
crpx_quasi_random_next_halton_original (crpx_quasi_random_t q)
{
  size_t i, j, k=0, iter_skip;
  q->iteration++;
  iter_skip = q->iteration;

  for (k = q->size - 1, i = 0; i < q->leap; i++) { // only enter this loop if q->size > HALTON_MAX_DIMENSION
    if (i < HALTON_MAX_DIMENSION + 1) iter_skip = (q->iteration) * halton_prime_numbers[i] + i; // "x prime" means "leap" and "+i" means "warm up"
    else iter_skip = (q->iteration) * i; // we're too far the original assumption of 1229 dimensions
    for (j = 1; j <= HALTON_MAX_DIMENSION; j++) {
      q->r[k--] = original_vdcorput (iter_skip, halton_prime_numbers[j]); // original GSL had prime_number[j] instead of j+1 but we use prime_number[0] to store "1"
    }
  }

  for (j = 1; j <= q->rem; j++) q->r[k--] = original_vdcorput (iter_skip, halton_prime_numbers[j]); // original GSL has prime_number[j] instead of j+1
  assert ((k == 0) && "not all quasi-random elements were generated in crpx_quasi_random_next_halton_original");
  return;
}

/* auxiliary functions */

static double
reverse_vdcorput (int x, int b)
{
  double r = 0., v = 1., binv = 1. / (double) b;
  while (x > 0) {
    v *= binv;
    r += v * (double) ((x % b == 0) ? 0 : b - (x % b));
    x /= b;
  }
  return r;
}

static double
original_vdcorput (int x, int b) // higher correlation, vdcorput is better
{
  double r = 0., v = 1., binv = 1. / (double) b;
  while (x > 0) {
    v *= binv;
    r += v * (double) (x % b);
    x /= b;
  }
  return r;
}

