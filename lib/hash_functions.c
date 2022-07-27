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


/*! \file hash_functions.c 
 *  \brief simple hash functions and random number generators from integers.
 */

#include "hash_functions.h"

void
crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4])
{
  FILE *fp;
  uint64_t x1[4]; // unsigned long int
  unsigned long long int x = 0; // courtesy of 32bits machines, they're not the same
  unsigned int i, idx, success;
  fp = fopen ("/dev/urandom", "r");
  if (fp == NULL) crpx_logger_info (cglob, "Could not open /dev/urandom to include it in the random seed\n");
  else {
    if (fread (seed, sizeof (uint64_t), 4, fp) != 1) crpx_logger_info (cglob, "Could not read from /dev/urandom to include it in the random seed\n");
    fclose (fp);
  }
  crpx_get_time_128bits (x1);
  seed[0] ^= x1[0];
  seed[1] ^= x1[1];
  seed[2] ^= x1[0];
  seed[3] ^= x1[1];
#ifdef HAVE_RDRND
  for (i=64, idx=0; idx < 4;i--) { // DRNG suggests 10 cycles 
    success = __builtin_ia32_rdrand64_step (&x);
    x1[idx] = x;
    idx += (success > 0); // avoid mispredicted branches (not performance-critical here though) 
  }
#endif
  x = (uint64_t) getpid() + getppid(); // process ID and parent process ID (can be small numbers)
  for (i=0; i<4;i++) seed[i] ^= (x1[i] + ((uint64_t)(i + 1) * x));

  return;
}

void
crpx_get_time_128bits (uint64_t time[2])
{
#if _POSIX_TIMERS
  struct timespec now;
  clock_gettime (CLOCK_REALTIME, &now);
  time[1] = now.tv_nsec; // always less than 1billion thus 32bits is enough
#else
  struct timeval now;
  gettimeofday (&now, NULL);
  time[1] = now.tv_usec; // always less than 1million thus 32bits is enough
#endif
  time[0] = now.tv_sec;
  return;
}

uint64_t 
crpx_wyhash64 (uint64_t *seed) // changes seed state (thus a PRNG) 
{  // https://github.com/lemire/testingRNG/blob/master/wyhash.c
  *seed += UINT64_C(0x60bee2bee120fc15);
  __uint128_t tmp;
  tmp = (__uint128_t)*seed * UINT64_C(0xa3b195354a39b70d);
  uint64_t m1 = (tmp >> 64) ^ tmp;
  tmp = (__uint128_t)m1 * UINT64_C(0x1b03738712fad5c9);
  uint64_t m2 = (tmp >> 64) ^ tmp;
  return m2;
}

uint64_t
crpx_splitmix64 (uint64_t *seed) // changes seed state (thus a PRNG) 
{ // https://github.com/lemire/testingRNG/blob/master/splitmix64.c
  uint64_t z = (*seed += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

uint64_t
crpx_fmix64 (uint64_t k) 
{ // https://github.com/lemire/clhash/blob/master/clhash.c ; mixer for murmurhash
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53ULL;
  k ^= k >> 33;
  return k;
}
