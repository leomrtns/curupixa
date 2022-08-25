/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html). 
 */

/*! \file random_number.c
 *  \brief high-level RNG, relying on lower level random_number_generators.c . */

#include "random_number.h"

inline double
biomcmc_rng_snorm (void)
{
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  if (biomcmc_random_number->have_rnorm64) {
    biomcmc_random_number->have_rnorm64 = false;
    return biomcmc_random_number->rnorm64;
  } 
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    /*    u and v are between -1 and +1 => -1 + 2 * U       */
    u = -1. + 2. * (biomcmc_rng_get_52 () / 4503599627370495.0); /* (2^52) - 1 => 1 included */
    v = -1. + 2. * (biomcmc_rng_get_52 () / 4503599627370495.0); /* (2^52) - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  biomcmc_random_number->rnorm64 = u * s;
  biomcmc_random_number->have_rnorm64 = true;
  return v * s;
}

/* in other places they use 
 * (genrand64_int64(context) >> 11) * (1.0/9007199254740991.0);         // [0,1]
 * (genrand64_int64(context) >> 11) * (1.0/9007199254740992.0);         // [0,1)
 * ((genrand64_int64(context) >> 12) + 0.5) * (1.0/4503599627370496.0); // (0,1)
 * in https://lemire.me/blog/2017/02/28/how-many-floating-point-numbers-are-in-the-interval-01/ he mentions 
 *    mantissa = 52 (23 for 32bits)  -->  random(0,2^53)/2^53 = [0,1) --> rnd() >> 11 / 9007199254740992
 *    and that rnd(0,1) * N is uniform(0,N) only for N < 2^53 
 */   

inline double
biomcmc_rng_unif (void)
{
  return (biomcmc_rng_get_52 () / 4503599627370495.0); /* 2^52 - 1 => 1 included */
}

inline double
biomcmc_rng_unif_pos (void)
{
  double x;
  do { x = (biomcmc_rng_get_52 () / 4503599627370495.0); } while (x < 2 * DBL_MIN);
  return x;
}

inline uint64_t
biomcmc_rng_unif_int64 (uint64_t n)
{
  uint64_t scale;
  uint64_t k;
  if (!n) biomcmc_error ("n must be larger than zero in uniform random number generator [64bits]");
  scale = 0xffffffffffffffffULL/n;

  do { k = biomcmc_rng_get ()/scale; } while (k >= n);
  return k;
}

inline double
biomcmc_rng_get_52 (void)
{
  /* In Matsumoto's MT19937 code they use 53 bits (total double precision) but I think that the integer fraction of a
   * double is only 52 - so the integer-to-double conversion should use only the first 52 bits */
  return (double) (biomcmc_rng_get () >> 12);
}


// return 0xffffffff & (h - (h >> 32)); // Fermat residue (https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c)
