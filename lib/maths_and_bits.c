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
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULLAR PURPOSE.  See the GNU General Public License for more 
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

inline uint32_t 
crpx_reverse_bits32 (uint32_t v) 
{ // https://graphics.stanford.edu/~seander/bithacks.html
  v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);// swap odd and even bits
  v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);// swap consecutive pairs
  v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);// swap nibbles ... 
  v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);// swap bytes
  v = ( v >> 16             ) | ( v               << 16);// swap 2-byte long pairs
  return v;
}

inline uint64_t 
crpx_interleave_64bits (uint64_t xylo) /* [y4y3y2y1 x4x3x2x1] --> [y4x4 y3x3 y2x2 y1x1] */
{ // https://github.com/yinqiwen/geohash-int/blob/master/geohash.c; https://graphics.stanford.edu/~seander/bithacks.html has for 32bits
  static const uint64_t B[] = { 0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF };
  static const unsigned int S[] = { 1, 2, 4, 8, 16 };
  uint64_t x = xylo & UINT64_C(0xFFFFFFFFFFFFFFFF);  // Interleave lower  bits of x and y, so the bits of x
  uint64_t y = xylo << 32;// are in the even positions and bits from y in the odd; //https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN

  // x and y must initially be less than 2**32. (@leomrtns: original algo uses two uint32_t xlo and ylo)
  x = (x | (x << S[4])) & B[4]; y = (y | (y << S[4])) & B[4];
  x = (x | (x << S[3])) & B[3]; y = (y | (y << S[3])) & B[3];
  x = (x | (x << S[2])) & B[2]; y = (y | (y << S[2])) & B[2];
  x = (x | (x << S[1])) & B[1]; y = (y | (y << S[1])) & B[1];
  x = (x | (x << S[0])) & B[0]; y = (y | (y << S[0])) & B[0];
  return x | (y << 1);
}

inline uint64_t 
crpx_deinterleave_64bits (uint64_t interleaved) /* [y4y3y2y1 x4x3x2x1] <-- [y4x4 y3x3 y2x2 y1x1] */
{ // https://github.com/yinqiwen/geohash-int/blob/master/geohash.c
  static const uint64_t B[] = { 0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF, 0x00000000FFFFFFFF };
  static const unsigned int S[] = { 0, 1, 2, 4, 8, 16 };
  uint64_t x = interleaved; ///reverse the interleave process (http://stackoverflow.com/questions/4909263/how-to-efficiently-de-interleave-bits-inverse-morton)
  uint64_t y = interleaved >> 1;
  x = (x | (x >> S[0])) & B[0]; y = (y | (y >> S[0])) & B[0];
  x = (x | (x >> S[1])) & B[1]; y = (y | (y >> S[1])) & B[1];
  x = (x | (x >> S[2])) & B[2]; y = (y | (y >> S[2])) & B[2];
  x = (x | (x >> S[3])) & B[3]; y = (y | (y >> S[3])) & B[3];
  x = (x | (x >> S[4])) & B[4]; y = (y | (y >> S[4])) & B[4];
  x = (x | (x >> S[5])) & B[5]; y = (y | (y >> S[5])) & B[5];
  return x | (y << 32);
}

/* https://github.com/lemire/Code-used-on-Daniel-Lemire-s-blog/blob/master/extra/search/shotgun/shotguntest.c binary search algos */

/* trick to swap a bit range within an integer, where r=swapped int; b=original int; n=length; i,j = positions to swap; x=temp var:
 *  x = ((b >> i) ^ (b >> j)) & ((1U << n) - 1); r = b ^ ((x << i) | (x << j));
 *  swaps the n consecutive bits starting at positions i and j (from the right) 
 *  See https://github.com/lemire/Code-used-on-Daniel-Lemire-s-blog/blob/master/extra/binomialcoef/binom.c for alternative */

inline int 
crpx_choose_n_k (int n, int k)
{
  if ((k > n) || (n > 100000)) return 0; // avoid overflow since 100000 choose 2 has more than 32 bits
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  /* curiosity: 512 choose 4 has almost 32 bits (thus should be a limit in quartets for instance) */
  int result = n;
  for (int i = 2; i <= k; ++i) { result *= (n-i+1); result /= i; }
  return result;
}

/* order-th lexicographically ordered set of k elements out of n; output is in result (of length at least k) */
void 
crpx_ordered_combination_n_k (int* result, int n, int k, int order)
{ // https://stackoverflow.com/questions/561/how-to-use-combinations-of-sets-as-test-data#794
  int r, j = 0, x = order+1; // order starts at zero, but lexico order starts at 1
  for (int i = 0; i < k-1; i++) {
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

