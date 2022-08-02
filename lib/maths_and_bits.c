/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be usefulull, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICullAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file maths_and_bits.c
 *  \brief Maths and bits functions.
 */

#include "maths_and_bits.h"

inline uint32_t
crpx_prev_power_of_two (uint32_t x)
{
  x |= (x >> 1); x |= (x >> 2); x |= (x >> 4); x |= (x >> 8); x |= (x >> 16);
  return x - (x >> 1);
}

inline uint32_t
crpx_next_power_of_two (uint32_t x)
{
  x--;
  x |= (x >> 1); x |= (x >> 2); x |= (x >> 4); x |= (x >> 8); x |= (x >> 16);
  return x + 1;
}

/* trick to swap a bit range within an integer, where r=swapped int; b=original int; n=length; i,j = positions to swap; x=temp var:
 *  x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); r = b ^ ((x << i) | (x << j));
 *  swaps the n consecutive bits starting at positions i and j (from the right) */

inline uint32_t 
crpx_choose_n_k (uint32_t n, uint32_t k)
{
  if ((k > n) || (n > 100000)) return 0; // avoid overflow since 100000 choose 2 has more than 32 bits
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  /* curiosity: 512 choose 4 has almost 32 bits (thus should be a limit in quartets for instance) */
  uint32_t result = n;
  for (uint32_t i = 2; i <= k; ++i) { result *= (n-i+1); result /= i; }
  return result;
}

/* order-th lexicographically ordered set of k elements out of n; output is in result (of length at least k) */
void 
crpx_ordered_combination_n_k (uint32_t* result, uint32_t n, uint32_t k, uint32_t order)
{ // https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
  uint32_t r, j = 0, x = order+1; // order starts at zero, but lexico order starts at 1
  for (uint32_t i = 0; i < k-1; i++) {
    result[i] = (i != 0) ? result[i-1] : 0;
    do {
      result[i]++;
      r = crpx_choose_n_k (n-result[i], k-(i+1));
      j += r;
    } while(j < x);
    j -= r;
  }
  result[k-1] = result[k-2] + x - j;
}

