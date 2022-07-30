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

size_t
crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size)
{
  size_t i = 0, j = 0, last = seed_size, first = 0;
  int success = 0; // gententropy() and syscall() return 0 on success and negative in failure (syscall returns errno)
  //uint16_t *seed_16bits = (uint16_t *)seed; // up to __builtin_ia32_rdrand64_step() but it tends to have too many zeroes
  uint32_t *seed_32bits = (uint32_t *)seed;

#ifdef HAVE_RDRND 
  for (i=0, j=0; (i < 2 * seed_size) && (j < (seed_size >> 2)); i++) { // DRNG suggests 10 attempts per rng, we do 2
//    success = __builtin_ia32_rdrand16_step (seed_16bits + j); // 16 bits = 2 bytes thus "j+=2" "j<size-1"
    success = __builtin_ia32_rdrand32_step (seed_32bits + j); // 64 bits = 8 bytes thus "j+=8" "j<size-1"
    j += (success > 0); // avoid mispredicted branches 
  } // it may ends with fewer than seed_size since not always succeed (and seed%2 may be > 0)
  j <<= 2;
  crpx_logger_verbose (cglob, "Random seeds produced by CPU crystal entropy (RDRAND): %lu", j);
#endif
  first = j;
  seed_size -= j;
  seed += j;
  for (i=0; (seed_size > 0) && (i < 2); i += (success<0)) { // i increments at every failure (which is a negative success)
    j = seed_size > 256 ? 256 : seed_size; /* maximum buffer size is 256 bytes */
#if (__GLIBC__ > 2 || __GLIBC_MINOR__ > 24)
    success = getentropy (seed, j);
#else
    success = syscall (SYS_getrandom, seed, j, 0);
#endif
    seed += j * (~success & 1);      // avoid mispredicted branches (not performance-critical here though)
    seed_size -= j * (~success & 1); // (success==0) is a success 
  }
  if (last-first) crpx_logger_verbose (cglob, "Random seeds produced by linux random: %lu", last - first);
  return last - seed_size;
}

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

#ifdef _POSIX_TIMERS
#define TIMEWARP 1.e9
#else
#define TIMEWARP 1.e6
#endif

double
crpx_update_elapsed_time_128bits (uint64_t past[2])
{
  uint64_t now[2];
  double seconds;
  crpx_get_time_128bits (now);
  if (now[1] < past[1]) seconds = (((double)(past[1] - now[1]) / (double)(TIMEWARP)) - 1. + (double)(now[0] - past[0]));
  else                  seconds = (((double)(now[1] - past[1]) / (double)(TIMEWARP))      + (double)(now[0] - past[0]));
  past[0] = now[0]; past[1] = now[1];
  return seconds;
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

uint64_t
crpx_hash_pearson (void *vkey, size_t len, const void *vseed) // seed must have >= 256 bytes
{ // https://github.com/maciejczyzewski/retter/blob/master/algorithms/Pearson/pearson.c
  uint8_t *key = (uint8_t *) vkey, *seed = (uint8_t *) vseed;
  uint64_t hash = 0;
  unsigned char h = 0;
  size_t i,j;
  for (j = 0; j < 8; j++) {
    h = seed[(key[0] + j) & 0xff ]; // 0xff = 255
    for (i = 1; i < len; i++) h = seed[h ^ key[i]];
    hash ^= h << (j * 8);
  }
  return hash;
}

uint32_t
crpx_hash_pseudocrc32 (uint32_t crc, void *vkey, size_t len, const uint32_t *seed) // seed must have >= 256 elements (of 32bits)
{ // modified from https://github.com/maciejczyzewski/retter/blob/master/algorithms/CRC/crc32.c (NOT CRC32)
  uint8_t *key = (uint8_t *) vkey;
  crc = crc ^ ~0U;
  while (len--) crc = seed[(crc ^ *key++) & 0xFF] ^ (crc >> 8);
  return crc ^ ~0U;
}

uint32_t 
crpx_hash_fletcher32 (uint16_t const *data, size_t words)
{ // https://github.com/maciejczyzewski/retter/tree/master/algorithms/Fletcher
/* The Fletcher checksum cannot distinguish between blocks of all 0  bits and blocks of all 1 bits. 
 * For example, if a 16-bit block in the data word changes from 0x0000 to 0xFFFF, the Fletcher-32 checksum remains the same. */
 uint32_t sum1 = 0xffff, sum2 = 0xffff;

 while (words) {
   unsigned tlen = words > 359 ? 359 : words;
   words -= tlen;
   do { sum2 += sum1 += *data++; } while (--tlen);
   sum1 = (sum1 & 0xffff) + (sum1 >> 16);
   sum2 = (sum2 & 0xffff) + (sum2 >> 16);
 }
 /* Second reduction step to reduce sums to 16 bits */
 sum1 = (sum1 & 0xffff) + (sum1 >> 16);
 sum2 = (sum2 & 0xffff) + (sum2 >> 16);
 return sum2 << 16 | sum1;
}

uint32_t
crpx_hash_jenkins (void *vkey, size_t len)
{ // https://github.com/maciejczyzewski/retter/tree/master/algorithms/Jenkins
  uint8_t *key = (uint8_t *) vkey;
  uint32_t hash, i;
  for(hash = i = 0; i < len; ++i) {
    hash += key[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
}
