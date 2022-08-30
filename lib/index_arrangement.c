
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

/*! \file index_arrangement.h 
 *  \brief combinations and permutations of indexes: a combination is just the subset of indices k out of n, and a permutation is each of the possible orderings of the subset.

 * based on the GSL (Copyright 1996~2007 Brian Gough under a GPL3) */ 

#include "index_arrangement.h"

/* permutations (size=3) : {0 1 2} {0 2 1} {1 0 2} {1 2 0} {2 0 1} {2 1 0} */

crpx_index_permutation_t 
crpx_index_permutation_new (crpx_global_t cglob, size_t n)
{
  crpx_index_permutation_t p = (crpx_index_permutation_t) crpx_malloc (cglob, sizeof (crpx_index_permutation_struct));
  if (!p) return NULL;
  p->size = n;
  p->idx = NULL;
  p->idx = (size_t *) crpx_malloc (cglob, n * sizeof (size_t));
  if (!p->idx) { crpx_free (cglob, p); return NULL; }
  crpx_link_add_global_pointer (cglob, p->cglob); // thread-safe increase of ref_counter
  crpx_index_permutation_reset (p);
  return p;
}

void
del_crpx_index_permutation (crpx_index_permutation_t p)
{
  if (!p) return;
  crpx_free (p->cglob, p->idx);  // free() with message if p->idx is NULL
  crpx_global_finalise (p->cglob); // it just decreases cglob->ref_counter unless order was inverted in main
  free (p);
}

void
crpx_index_permutation_reset (crpx_index_permutation_t p)
{
  if (!p) return;
  for (size_t i = 0; i < p->size; i++) p->idx[i] = i;
}

bool
crpx_index_permutation_swap (crpx_index_permutation_t p, const size_t i, const size_t j)
{
  if (i >= p->size) {
    crpx_logger_error (p->cglob, "index_permutation_swap: location i = %zu out of range (larger than size %zu)", i, p->size);
    return false;
  }
  if (j >= p->size) {
    crpx_logger_error (p->cglob, "index_permutation_swap: location j = %zu out of range (larger than size %zu)", j, p->size);
    return false;
  }

  if (i != j) {
    size_t tmp = p->idx[i];
    p->idx[i] = p->idx[j];
    p->idx[j] = tmp;
  }
  return true;
}

void
crpx_index_permutation_swap_no_checks (crpx_index_permutation_t p, const size_t i, const size_t j)
{
  size_t tmp = p->idx[i];
  p->idx[i] = p->idx[j];
  p->idx[j] = tmp;
}

bool
crpx_index_permutation_is_valid (const crpx_index_permutation_t p) // just a check, therefore doesn't raise an error but logs a debug
{
  size_t i, j;
  for (i = 0; i < p->size; i++) {
    if (p->idx[i] >= p->size) { 
      crpx_logger_debug (p->cglob, "index_permutation check: index %zu at location %zu is out of range (larger than size %zu)", p->idx[i], i, p->size);
      return false;
    }
    for (j = 0; j < i; j++) if (p->idx[i] == p->idx[j]) {
      crpx_logger_debug (p->cglob, "index_permutation check: index %zu is repeated at locations %zu and %zu", p->idx[i], i, j);
      return false;
    }
  }
  return true;
}

void
crpx_index_permutation_reverse (crpx_index_permutation_t p) // reverse all elements backwards
{
  size_t i, j, tmp, halfsize = p->size / 2;
  for (i = 0; i < halfsize; i++) {
    j = p->size - i - 1;
    tmp = p->idx[i] ;
    p->idx[i] = p->idx[j] ;
    p->idx[j] = tmp;
  }
}

bool 
crpx_index_permutation_inverse (crpx_index_permutation_t q, const crpx_index_permutation_t p)
{ // q[] will have the elements of p[] which would make p[] the identity index (i.e. p[q[i]] = i)
  size_t i ;
  if (q->size != p->size) {
    crpx_logger_error (p->cglob, "index_permutation_inverse: q->size = %zu, p->size = %zu", q->size, p->size);
    return false;
  }
  for (i = 0; i < q->size; i++) q->idx[p->idx[i]] = i;
  return true;
}

bool
crpx_index_permutation_next (crpx_index_permutation_t p)
{ // Replaces p with the next permutation (in lexico ordering). Returns false if there is no next permutation (i.e. iterated over all permutations).
  size_t i, j, k, tmp;

  if (p->size < 2) return false;
  i = p->size - 2;
  while ((p->idx[i] > p->idx[i+1]) && (i)) i--;
  if ((i == 0) && (p->idx[0] > p->idx[1])) return false;
  k = i + 1;
  for (j = i + 2; j < p->size; j++ ) { if ((p->idx[j] > p->idx[i]) && (p->idx[j] < p->idx[k])) k = j; }

  tmp = p->idx[i]; /* swap i and k */
  p->idx[i] = p->idx[k];
  p->idx[k] = tmp;

  for (j = i + 1; j <= ((p->size + i) / 2); j++) {
    tmp = p->idx[j];
    p->idx[j] = p->idx[p->size + i - j];
    p->idx[p->size + i - j] = tmp;
  }
  return true;
}

bool
crpx_index_permutation_prev (crpx_index_permutation_t p)
{
  size_t i, j, k, tmp;

  if (p->size < 2) return false;
  i = p->size - 2;
  while ((p->idx[i] < p->idx[i+1]) && (i)) i--;
  if ((i == 0) && (p->idx[0] < p->idx[1])) return false;
  k = i + 1;
  for (j = i + 2; j < p->size; j++ ) { if ((p->idx[j] < p->idx[i]) && (p->idx[j] > p->idx[k])) k = j; }

  tmp = p->idx[i];  /* swap i and k */
  p->idx[i] = p->idx[k];
  p->idx[k] = tmp;

  for (j = i + 1; j <= ((p->size + i) / 2); j++) {
    tmp = p->idx[j];
    p->idx[j] = p->idx[p->size + i - j];
    p->idx[p->size + i - j] = tmp;
  }
  return true;
}

bool
crpx_index_permutation_combine (crpx_index_permutation_t p, const crpx_index_permutation_t pa, const crpx_index_permutation_t pb)
{
  size_t i;
  if (pa->size != p->size) {
    crpx_logger_error (p->cglob, "index_permutation_combine: input pa->size = %zu but output size = %zu", pa->size, p->size);
    return false;
  }
  if (pb->size != p->size) {
    crpx_logger_error (p->cglob, "index_permutation_combine: input pb->size = %zu but output size = %zu", pb->size, p->size);
    return false;
  }
  for (i = 0; i < p->size; i++) p->idx[i] = pb->idx[ pa->idx[i] ];
  return true;
}

/* combination (k = 2, n = 4) : { 0 1 }{ 0 2 }{ 0 3 }{ 1 2 }{ 1 3 }{ 2 3 } notice that all are in increasing order */ 

crpx_index_combination_t 
crpx_index_combination_new (crpx_global_t cglob, size_t n, size_t k)
{
  crpx_index_combination_t c = (crpx_index_combination_t) crpx_malloc (cglob, sizeof (crpx_index_combination_struct));
  if (!c) return NULL;
  c->n = n;
  c->k = k;
  c->idx = NULL;
  c->idx = (size_t *) crpx_malloc (cglob, k * sizeof (size_t));
  if (!c->idx) { crpx_free (cglob, c); return NULL; }
  crpx_link_add_global_pointer (cglob, c->cglob); // thread-safe increase of ref_counter
  crpx_index_combination_reset_first (c);
  return c;
}

void
del_crpx_index_combination (crpx_index_combination_t c)
{
  if (!c) return;
  crpx_free (c->cglob, c->idx);
  crpx_global_finalise (c->cglob); // it just decreases cglob->ref_counter unless order was inverted in main
  free (c);
}

void
crpx_index_combination_reset_first (crpx_index_combination_t c)
{
  if (!c) return;
  for (size_t i = 0; i < c->k; i++) c->idx[i] = i;
}

void
crpx_index_combination_reset_last (crpx_index_combination_t c)
{
  if (!c) return;
  for (size_t i = 0; i < c->k; i++) c->idx[i] = c->n - c->k + i;
}

bool 
crpx_index_combination_is_valid (crpx_index_combination_t c)
{
  size_t i, j, ci;
  if(c->k > c->n) {
    crpx_logger_debug (c->cglob, "index_combination check: k=%zu is larger than size %zu)", c->k, c->n);
    return false;
  }
  for (i = 0; i < c->k; i++)  {
    ci = c->idx[i];
    if (ci >= c->n) {
      crpx_logger_debug (c->cglob, "index_combination check: index %zu at location %zu outside range (larger than size %zu)", ci, i, c->n);
      return false;
    }
    for (j = 0; j < i; j++) {
      if (c->idx[j] == ci) {
        crpx_logger_debug (c->cglob, "index_combination check: index %zu is duplicate, at locations %zu and %zu", ci, i, j);
        return false;
      }
      if (c->idx[j] > ci) {
        crpx_logger_debug (c->cglob, "index_combination check: idx[%zu]=%zu > idx[%zu]=%zu not in increasing order ", j, c->idx[j], i, ci);
        return false;
      }
    }
  }
  return true; 
}

bool
crpx_index_combination_next (crpx_index_combination_t c)
{ // Replaces c with the next combination (in lexico ordering). Returns false if all combinations were given
  if (c->k == 0) return false;
  size_t i = c->k - 1;

  while ((i > 0) && (c->idx[i] == (c->n - c->k + i))) i--;
  if ((i == 0) && (c->idx[i] == (c->n - c->k))) return false;
  c->idx[i]++;
  for (; i < c->k - 1; i++) c->idx[i + 1] = c->idx[i] + 1;
  return true;
}

bool
crpx_index_combination_prev (crpx_index_combination_t c)
{ // Replaces c with the previous combination (in lexico ordering). Returns false if all combinations were given
  if (c->k == 0) return false;
  size_t i = c->k - 1;

  while ((i > 0) && (c->idx[i] == (c->idx[i-1] + 1))) i--;
  if ((i == 0) && (c->idx[i] == 0)) return false;
  c->idx[i++]--;
  for(; i < c->k; i++) c->idx[i] = c->n - c->k + i;
  return true;
}

