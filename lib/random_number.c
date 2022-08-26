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

void
crpx_set_random_generator (crpx_global_t cglob, uint8_t rng_id, uint64_t seed)
{
  char rname[32] = {'\0'};
  switch (rng_id) {
    case 0:  cglob->rng_get = &crpx_rng_wyhash_state64;       cglob->rng_size = 1; strcpy (rname, "wyhash_64"); break;
    case 1:  cglob->rng_get = &crpx_rng_lehmer_seed128;       cglob->rng_size = 2; strcpy (rname, "lehmer_64"); break;
    case 2:  cglob->rng_get = &crpx_rng_moremur_seed64;       cglob->rng_size = 1; strcpy (rname, "moremur_64"); break;
    case 3:  cglob->rng_get = &crpx_rng_splitmix_seed64;      cglob->rng_size = 1; strcpy (rname, "splitmix_64"); break;
    case 4:  cglob->rng_get = &crpx_rng_romu_seed128;         cglob->rng_size = 2; strcpy (rname, "romu_128"); break;
    case 5:  cglob->rng_get = &crpx_rng_jenkins13_seed256;    cglob->rng_size = 4; strcpy (rname, "jenkins13_256"); break;
    case 6:  cglob->rng_get = &crpx_rng_jenkins19_seed256;    cglob->rng_size = 4; strcpy (rname, "jenkins19_256"); break;
    case 7:  cglob->rng_get = &crpx_rng_xorshift_star_seed64; cglob->rng_size = 1; strcpy (rname, "xorshift_s_64"); break;
    case 8:  cglob->rng_get = &crpx_rng_romu_seed192;         cglob->rng_size = 3; strcpy (rname, "romu_192"); break;
    case 9:  cglob->rng_get = &crpx_xoroshiro_pv6_seed128;    cglob->rng_size = 2; strcpy (rname, "xoro_pv6_128"); break;
    case 10: cglob->rng_get = &crpx_xoroshiro_pv8_seed128;    cglob->rng_size = 2; strcpy (rname, "xoro_pv8_128"); break;
    case 11: cglob->rng_get = &crpx_rng_romu_seed256;         cglob->rng_size = 4; strcpy (rname, "romu_256"); break;
    case 12: cglob->rng_get = &crpx_rng_xorshift_p_seed128;   cglob->rng_size = 2; strcpy (rname, "xorshift_p_128"); break;
    case 13: cglob->rng_get = &crpx_rng_rrmixer_seed64;       cglob->rng_size = 1; strcpy (rname, "rrmixer_64"); break;
    case 14: cglob->rng_get = &crpx_xoroshiro_pp_seed128;     cglob->rng_size = 2; strcpy (rname, "xoro_pp_128"); break;
    case 15: cglob->rng_get = &crpx_xoroshiro_star_seed256;   cglob->rng_size = 4; strcpy (rname, "xoro_s_256"); break;
    case 16: cglob->rng_get = &crpx_rng_wyrand_seed64;        cglob->rng_size = 1; strcpy (rname, "wyrand_64"); break;
    case 17: cglob->rng_get = &crpx_xoroshiro_pp_seed256;     cglob->rng_size = 4; strcpy (rname, "xoro_pp_256"); break;
    case 18: cglob->rng_get = &crpx_rng_pcg_seed256;          cglob->rng_size = 4; strcpy (rname, "pcg_256"); break;
    default: cglob->rng_get = &crpx_rng_mt19937_seed2504;   cglob->rng_size = 313; strcpy (rname, "mt19937"); break;
  }
  
  size_t n_bytes = cglob->rng_size * cglob->nthreads * sizeof (uint64_t), success_bytes = 0;
  cglob->rng_seed_vector = (uint64_t *) crpx_realloc (n_bytes); // malloc() if first time; we neglect current contents of the vector
  uint8_t *seed_vector_bytes = (uint8_t *) cglob->rng_seed_vector;
  if (!seed) success_bytes = crpx_generate_bytesized_random_seeds_from_cpu (cglob, cglob->rng_seed_vector, n_bytes);
  if (success_bytes < n_bytes)  // or seed==0 or cpu random was not successful
    crpx_generate_bytesized_random_seeds_from_seed (cglob, seed_vector_bytes + success_bytes, n_bytes - success_bytes, seed);

  crpx_logger_verbose (cglob, "Random number generator set to '%s' (using %u bytes of state)", rname, cglob->rng_size * sizeof (uint64_t));
}

inline uint64_t
crpx_random_64bits (crpx_global_t cglob)
{
  return cglob->rng_get(cglob->rng_seed_vector + cglob->rng_size * CRPX_THREAD_NUM);
}

inline uint32_t
crpx_random_32bits (crpx_global_t cglob)
{
  uint64_t h = crpx_random_64bits (cglob);
  return (uint32_t)(0xffffffff & (h - (h >> 32))); // Fermat residue (https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c)
}

inline uint64_t
crpx_random_range (crpx_global_t cglob, uint64_t n)
{
  uint64_t scale, k;
  assert (n > 0 && "n must be positive otherwise will lead to division by zero");
  if (!n) {
    n = 1;
    crpx_logger_warning (cglob, "curupixa_random_range(n=0) is not defined, will assume n=1");
  }
  scale = UINT64_MAX / n;
  do { k = crpx_random_64bits(cglob)/scale; } while (k >= n);
  return k;
}

/* in https://lemire.me/blog/2017/02/28/how-many-floating-point-numbers-are-in-the-interval-01/ he mentions 
 *    mantissa = 52 (23 for 32bits)  -->  random(0,2^53)/2^53 = [0,1) --> rnd() >> 11 / 9007199254740992
 *    and that rnd(0,1) * N is uniform(0,N) only for N < 2^53   */   
inline double
crpx_random_double (crpx_global_t cglob) // [0,1)
{
  return (double)(crpx_random_64bits(cglob) >> 11) / 9007199254740992.0; // 2^53
}

inline double
crpx_random_double_include_one (crpx_global_t cglob) // [0,1]
{
  return (double)(crpx_random_64bits(cglob) >> 11) / 9007199254740991.0; // 2^53-1
}

inline double
crpx_random_double_positive (crpx_global_t cglob) // (0,1)
{
  double x;
  do { x = (crpx_random_64bits(cglob) >> 11) / 9007199254740992.0; } while (x < 2 * DBL_MIN);
 // return ((double)(genrand64_int64(context) >> 12) + 0.5) / 4503599627370496.0; // alternative 
}

inline double
crpx_random_double_positive_include_one (crpx_global_t cglob) // (0,1]
{
  double x;
  do { x = (crpx_random_64bits(cglob) >> 11) / 9007199254740991.0; } while (x < 2 * DBL_MIN);
}

inline double
crpx_random_normal (crpx_global_t cglob, double *extra_result) // generates _2_ random numbers, one is stored in extra_result 
{ // from biomcmc
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    /*    u and v are between -1 and +1 => -1 + 2 * U       */
    u = -1. + 2. * ((double)(crpx_random_64bits (cglob) >> 11) / 9007199254740991.0); /* (2^53) - 1 => 1 included */
    v = -1. + 2. * ((double)(crpx_random_64bits (cglob) >> 11) / 9007199254740991.0); /* (2^53) - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  (*extra_result) = u * s;
  return v * s;
}
