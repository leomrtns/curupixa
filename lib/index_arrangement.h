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
 *  \brief combinations C(M) and permutations P(k,N) for an index vector (i.e. assuming indices 0...M-1 or 0...k-1)
 *  Combinations are 0...M-1 in all orderings, and permutations are all k numbers out of N in increasing order */ 

#ifndef _curupixa_index_arrangement_h_
#define _curupixa_index_arrangement_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "random_number.h"

typedef struct {
  size_t size, *idx; 
  crpx_global_t cglob;
} crpx_index_permutation_struct, *crpx_index_permutation_t;

typedef struct {
  size_t k, n, *idx; 
  crpx_global_t cglob;
} crpx_index_combination_struct, *crpx_index_combination_t;

/* permutations (size=3) : {0 1 2} {0 2 1} {1 0 2} {1 2 0} {2 0 1} {2 1 0} */

crpx_index_permutation_t crpx_index_permutation_new (crpx_global_t cglob, size_t n);
void del_crpx_index_permutation (crpx_index_permutation_t p);
void crpx_index_permutation_reset (crpx_index_permutation_t p);
bool crpx_index_permutation_swap (crpx_index_permutation_t p, const size_t i, const size_t j);
void crpx_index_permutation_swap_no_checks (crpx_index_permutation_t p, const size_t i, const size_t j);
bool crpx_index_permutation_is_valid (const crpx_index_permutation_t p); // just a check, therefore doesn't raise an error but logs a debug
void crpx_index_permutation_reverse (crpx_index_permutation_t p); // reverse all elements backwards
bool crpx_index_permutation_inverse (crpx_index_permutation_t q, const crpx_index_permutation_t p);
bool crpx_index_permutation_next (crpx_index_permutation_t p); 
bool crpx_index_permutation_prev (crpx_index_permutation_t p);
bool crpx_index_permutation_combine (crpx_index_permutation_t p, const crpx_index_permutation_t pa, const crpx_index_permutation_t pb);

/* combination (k = 2, n = 4) : { 0 1 }{ 0 2 }{ 0 3 }{ 1 2 }{ 1 3 }{ 2 3 } notice that all are in increasing order */ 

crpx_index_combination_t crpx_index_combination_new (crpx_global_t cglob, size_t n, size_t k);
void del_crpx_index_combination (crpx_index_combination_t c);
void crpx_index_combination_reset_first (crpx_index_combination_t c);
void crpx_index_combination_reset_last (crpx_index_combination_t c);
bool crpx_index_combination_is_valid (crpx_index_combination_t c);
bool crpx_index_combination_next (crpx_index_combination_t c);
bool crpx_index_combination_prev (crpx_index_combination_t c);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
