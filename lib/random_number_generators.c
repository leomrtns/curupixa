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

/*! \file random_number_generators.c
 *  \brief 64 bits random number generators from integers. */

#include "random_number_generators.h"
#include "internal_random_constants.h" // not available to the user, only locally

uint64_t 
crpx_rng_wyhash_state64 (void *vstate)
{  // https://github.com/lemire/testingRNG/blob/master/wyhash.c
  uint64_t *state = (uint64_t *)vstate;
  *state += UINT64_C(0x60bee2bee120fc15);
  __uint128_t tmp;
  tmp = (__uint128_t)*state * UINT64_C(0xa3b195354a39b70d);
  uint64_t m1 = (tmp >> 64) ^ tmp;
  tmp = (__uint128_t)m1 * UINT64_C(0x1b03738712fad5c9);
  uint64_t m2 = (tmp >> 64) ^ tmp;
  return m2;
}

uint64_t
crpx_rng_splitmix_seed64 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/splitmix64.c
  uint64_t *state = (uint64_t *) vstate;
  uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

uint64_t 
crpx_lehmer_seed128 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
  __uint128_t *state = (__uint128_t *) vstate;
  *state *= UINT64_C(0xda942042e4dd58b5);
  return (uint64_t) (*state >> 64);
}

uint64_t 
crpx_wyrand_seed64 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/source/wyrand.h
  uint64_t *state = (uint64_t *) vstate;
  *state += UINT64_C(0xa0761d6478bd642f);
  __uint128_t t = (__uint128_t)(*state) * ((*state) ^ UINT64_C(0xe7037ed1a0b428db));
  return (uint64_t)((t >> 64) ^ t);
}

uint64_t 
crpx_rng_jenkins13_seed256 (void *vstate) // 13 bits of avalanche
{ // http://burtleburtle.net/bob/rand/smallprng.html see "64 bits variants" (and testingRNG below) 
  uint64_t *s = (uint64_t *) vstate; 
  uint64_t e = s[0] - ROTL64(s[1], 39);
  s[0] = s[1] ^ ROTL64(s[2], 11);
  s[1] = s[2] + s[3];
  s[2] = s[3] + e;
  s[3] = e + s[0];
  return s[3];
}

uint64_t 
crpx_rng_jenkins19_seed256 (void *vstate) // 18.4 bits of avalanche
{ // https://github.com/lemire/testingRNG/blob/master/source/jenkinssmall.h
  uint64_t *s = (uint64_t *) vstate; 
  uint64_t e = s[0] - ROTL64(s[1], 7);
  s[0] = s[1] ^ ROTL64(s[2], 13);
  s[1] = s[2] + ROTL64(s[3], 37);
  s[2] = s[3] + e;
  s[3] = e + s[0];
  return s[3];
}

uint64_t
crpx_rng_rrmixer_seed64 (void *vstate)
{ // general hash-to-rng trick: increment state and hash it 
  uint64_t *state = (uint64_t *) vstate;
  uint64_t k = (*state += UINT64_C(0x2adca2f5d6da1507)); // large prime by @leomrtns
  k ^= ROTR64(k, 49) ^ ROTR64(k, 24);
  k *= 0x9fb21c651e98df25LL;
  k ^= k >> 28;
  k *= 0x9fb21c651e98df25LL;
  return k ^ k >> 28;
}

uint64_t
crpx_rng_moremur_seed64 (void *vstate)
{ 
  uint64_t *state = (uint64_t *) vstate;
  uint64_t x = (*state += UINT64_C(0x0be40fe266ab1ec7)); // large prime by @leomrtns
  x ^= x >> 27;
  x *= 0x3C79AC492BA7B653UL;
  x ^= x >> 33;
  x *= 0x1C69B3F74AC4AE35UL;
  x ^= x >> 27;
  return x;
}

uint64_t
crpx_rng_romu_seed256 (void *vstate) // romu_quad: 4 x uint64_t 
{ // https://github.com/opencoff/portable-lib/blob/master/src/romu-rand.c 
  uint64_t *r = (uint64_t *) vstate;
  uint64_t x = r[0], y = r[1], z = r[2], w=r[3];
  r[3] = 15241094284759029579ULL * z;
  r[0] = z + ROTL64(w, 52);
  r[1] = y - x;
  r[2] = ROTL64(y+w, 19);
  return x;
}

uint64_t
crpx_rng_romu_seed192 (void *vstate) // romu_trio: 3 x uint64_t 
{ // https://github.com/opencoff/portable-lib/blob/master/src/romu-rand.c 
  uint64_t *r = (uint64_t *) vstate;
  uint64_t x = r[0], y = r[1], z = r[2];
  r[0] = 15241094284759029579ULL * z;
  r[1] = ROTL64(y - x, 12);
  r[2] = ROTL64(z - y, 44);
  return x;
}

uint64_t
crpx_rng_romu_seed128 (void *vstate) // romu_duo: 2 x uint64_t 
{ // https://github.com/opencoff/portable-lib/blob/master/src/romu-rand.c 
  uint64_t *r = (uint64_t *) vstate;
  uint64_t x = r[0];
  r[0] = 15241094284759029579u * r[1];
  r[1] = ROTL64(r[1], 36) + ROTL64(r[1], 15) - x;
  return x;
}

/* See also http://prng.di.unimi.it/xoroshiro128plusplus.c and http://xoroshiro.di.unimi.it/xoroshiro128plus.c 
 * https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c contains others, and jump functions */
uint64_t
crpx_xoroshiro_pv6_seed128 (void *vstate) // 128+ V 2016
{ // https://github.com/opencoff/portable-lib/blob/master/src/xoroshiro.c
  uint64_t *v = (uint64_t *) vstate;
  uint64_t v0 = v[0], v1 = v[1];
  uint64_t result = v[0] + v[1]; 
  v1 ^= v0; v[0] = ROTL64(v0, 55) ^ v1 ^ (v1 << 14); v[1] = ROTL64(v1, 36);
  return result;
}

uint64_t
crpx_xoroshiro_pv8_seed128 (void *vstate) // 128+ V 2018
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
  uint64_t v0 = v[0], v1 = v[1];
  uint64_t result = v[0] + v[1]; 
  v1 ^= v0; v[0] = ROTL64(v0, 24) ^ v1 ^ (v1 << 16); v[1] = ROTL64(v1, 37);
  return result;
}

uint64_t
crpx_xoroshiro_pp_seed128 (void *vstate) // 128++
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
  uint64_t v0 = v[0], v1 = v[1];
  uint64_t result = ROTL64(v0 + v1, 17) + v0; 
  v1 ^= v0; v[0] = ROTL64(v0, 49) ^ v1 ^ (v1 << 21); v[1] = ROTL64(v1, 28);
  return result;
}

uint64_t
crpx_xoroshiro_pp_seed256 (void *vstate) 
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
  uint64_t t = v[1] << 17; 
  uint64_t result = ROTL64(v[1] * 5, 7) * 9;
  v[2] ^= v[0]; v[3] ^= v[1]; v[1] ^= v[2]; v[0] ^= v[3]; v[2] ^= t; v[3] = ROTL64(v[3], 45);
  return result;
}

uint64_t
crpx_xoroshiro_star_seed256 (void *vstate) 
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
  uint64_t t = v[1] << 17; 
  uint64_t result = ROTL64(v[0] + v[3], 23) + v[0];
  v[2] ^= v[0]; v[3] ^= v[1]; v[1] ^= v[2]; v[0] ^= v[3]; v[2] ^= t; v[3] = ROTL64(v[3], 45);
  return result;
}

uint64_t
crpx_xorshift_star_seed64 (void *vstate)
{ // https://github.com/opencoff/portable-lib/blob/master/src/xorshift.c
  uint64_t *s = (uint64_t *) vstate;  
  *s ^= *s >> 12; *s ^= *s << 25; *s ^= *s << 27;
  return (*s) * 2685821657736338717ULL;
}

uint64_t
crpx_xorshift_p_seed128 (void *vstate)
{ // https://github.com/opencoff/portable-lib/blob/master/src/xorshift.c
  uint64_t *s = (uint64_t *) vstate;
  uint64_t v1 = s[0], v0 = s[1];
  s[0] = v0; 
  v1 ^= (v1 << 23);
  s[1] = v1 ^ v0 ^ (v1 >> 18) ^ (v0 >> 5);
  return s[1] + v0;
}

#define PCG_DEFAULT_MULTIPLIER_128 ((((__uint128_t)2549297995355413924ULL) << 64) + 4865540595714422341ULL)
uint64_t
crpx_pcg_seed256 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/source/pcg64.h
  __uint128_t *s = (__uint128_t *) vstate;
  uint64_t value;
  unsigned int rot;
  s[0] = s[0] * PCG_DEFAULT_MULTIPLIER_128 + s[1]; // s[1] never changes but _must_be_ odd;
  value = ((uint64_t)(s[0] >> 64u)) ^ (uint64_t) s[0];
  rot = s[0] >> 122u;
  return (value >> rot) | (value << ((-rot) & 63));
}

/* 32 bits */

uint32_t
crpx_rng_abyssinian_seed128 (void *vstate) // 2 x uint64_t 
{ // https://github.com/opencoff/portable-lib/blob/master/src/abyssinian_rand.c
  uint64_t *state = (uint64_t *)vstate;
  state[0] = (uint64_t)0xfffd21a7 * (uint32_t)state[0] + (uint32_t)(state[0] >> 32);
  state[1] = (uint64_t)0xfffd1361 * (uint32_t)state[1] + (uint32_t)(state[1] >> 32);
  return ROTL32((uint32_t)state[0], 7) + (uint32_t)state[1];
}

uint32_t 
crps_rng_widynski_seed192 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/source/widynski.h
  uint64_t *state = (uint64_t *)vstate;
  state[0] *= state[0];
  state[0] += (state[1] += state[2]);
  state[0] = (state[0] >> 32) | (state[0] << 32);
  return state[0];
}

uint32_t 
crpx_rng_jenkins8_seed128 (void *vstate) // 8 bits of avalanche 
{ // http://burtleburtle.net/bob/rand/smallprng.html amd testingRNG/source/jenkinssmall.h
  uint32_t *s = (uint32_t *) vstate; 
  // other combinations: (9,16), (9,24), (10,16), (10,24), (11,16), (11,24), (25,8), (25,16), (26,8), (26,16), (26,17), (27,16)
  uint32_t e = s[0] - ROTL32(s[1], 27);
  s[0] = s[1] ^ ROTL32(s[2], 17);
  s[1] = s[2] + s[3];
  s[2] = s[3] + e;
  s[3] = e + s[0];
  return s[3];
}

uint32_t 
crpx_rng_jenkins13_seed128 (void *vstate) // 13 bits of avalanche
{ // https://github.com/lemire/testingRNG/blob/master/source/jenkinssmall.h
  uint32_t *s = (uint32_t *) vstate; 
  // other combinations: (3,14,24), (3,25,15), (4,15,24), (6,16,28), (7,16,27), (8,14,3), (11,16,23), (12,16,22), (12,17,23), (13,16,22), (15,25,3), 
  //  (16,9,3), (17,9,3), (17,27,7), (19,7,3), (23,15,11), (23,16,11), (23,17,11), (24,3,16), (24,4,16), (25,14,3), (27,16,6), (27,16,7).
  uint32_t e = s[0] - ROTL32(s[1], 23);
  s[0] = s[1] ^ ROTL32(s[2], 16);
  s[1] = s[2] + ROTL32(s[3], 11);
  s[2] = s[3] + e;
  s[3] = e + s[0];
  return s[3];
}

uint64_t
crpx_rng_mt19937_seed2504 (void *state) // needs 312 uint64_t for random state and last one is a counter
{ // adapted from biomcmc
  static const uint64_t mag01[2]={ 0ULL, 0xB5026F5AA96619E9ULL}; /* this is magic vector, don't change */
  uint64_t *r = (uint64_t *) state;
  uint64_t x;

  if (r[312] >= 312) { /* generate all 312 words at once; r[312] is counter */
    int i;
    for (i = 0; i < 156; i++) {
      x = (r[i] & 0xFFFFFFFF80000000ULL)| (r[i+1] & 0x7FFFFFFFULL);
      r[i] = r[i+156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    for (; i < 311; i++) {
      x = (r[i] & 0xFFFFFFFF80000000ULL) | (r[i+1] & 0x7FFFFFFFULL);
      r[i] = r[i-156] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    }
    x = (r[311] & 0xFFFFFFFF80000000ULL) | (r[0] & 0x7FFFFFFFULL);
    r[311] = r[155] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
    r[312] = 0; // zero counter
  }

  x = r[ r[312]++ ];
  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);
  return x;
}


/* Generation of seeds */ 

void 
crpx_pcg_set_seed256 (void *vstate, uint64_t seed)
{ // https://github.com/lemire/testingRNG/blob/master/source/pcg64.h
  uint64_t *s = (uint64_t *) vstate;
  uint64_t s0 = seed; // will modify seed's state so we use a copy
  if (!s0) s0 = 0x1576359c7fcbd9dfULL; // arbitrary prime number by @leomrtns
  __uint128_t initstate, initseq;
  s[0] = crpx_rng_splitmix_seed64 (&s0); s0++;
  s[1] = crpx_rng_splitmix_seed64 (&s0); s0++;
  s[2] = crpx_rng_splitmix_seed64 (&s0); s0++;
  s[3] = crpx_rng_splitmix_seed64 (&s0); s0++;
  initstate = ((((__uint128_t)s[0]) << 64) + s[1]);
  initseq   = ((((__uint128_t)s[2]) << 64) + s[3]);
  initseq = (initseq << 1u) | 1u; // makes sure it's odd
  initstate += initseq;
  initstate = initstate * PCG_DEFAULT_MULTIPLIER_128 + initseq;
  s[0] = (uint64_t) (initstate >> 64); s[1] = (uint64_t) initstate;
  s[2] = (uint64_t) (initseq >> 64);   s[3] = (uint64_t) initseq;
}

void
cprx_rng_abyssinian_set_seed128 (void *vstate, uint32_t seed)
{ // https://github.com/opencoff/portable-lib/blob/master/src/abyssinian_rand.c
  uint64_t *state = (uint64_t *)vstate;
  if (!seed) seed = 0xc3fdc7fU; // arbitrary prime number by @leomrtns
  uint32_t seed_2 = seed;
  seed += seed_2;
  seed_2 += seed;
  uint64_t seed_x = 0x9368e53c2f6af274ULL ^ seed;
  uint64_t seed_y = 0x586dcd208f7cd3fdULL ^ seed_2;

  seed_x *= 0xff51afd7ed558ccdULL;// Based on the mixing functions of MurmurHash3
  seed_x ^= seed_x >> 33;
  seed_x *= 0xc4ceb9fe1a85ec53ULL;
  seed_x ^= seed_x >> 33;

  seed_y *= 0xff51afd7ed558ccdULL;
  seed_y ^= seed_y >> 33;
  seed_y *= 0xc4ceb9fe1a85ec53ULL;
  seed_y ^= seed_y >> 33;
  // Inlined Next(): Discard first output
  seed_x = (uint64_t)0xfffd21a7 * (uint32_t)seed_x + (uint32_t)(seed_x >> 32);
  seed_y = (uint64_t)0xfffd1361 * (uint32_t)seed_y + (uint32_t)(seed_y >> 32);
  state[0] = seed_x;
  state[1] = seed_y;
}

void
crpx_rng_mt19937_set_seed2504 (void *state, uint64_t seed)
{ // adapted from biomcmc
  uint64_t *r = (uint64_t *) state;
  uint64_t s0 = seed; 
  if (!s0) s0 = 0x8fc18365c966079ULL; // arbitrary prime number by @leomrtns
  r[312] = 313; // counter, forcing generation of 312 words at first call

  for (int i = 0; i < 312; i++) r[i] = crpx_rng_splitmix_seed64 (&s0);
}

