/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT AND Apache-2.0 AND BSD-3-Clause AND Unlicense */
/*! \file random_number_generators.c
    Small stateful generators, modified for explicit Curupixa state and C17. Original permissive sources and notices are recorded in docs/random-hashing.md.
    Blackman and Vigna's xoshiro256++ is the default; the other generators retain their historical identifiers. */
#include "random_number_generators.h"
#include "hash_functions_generators.h"
#include "internal_random_constants.h"

uint64_t 
crpx_rng_wyhash_state64 (void *vstate)
{
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
{
  uint64_t *state = (uint64_t *) vstate;
  uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
  z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
  return z ^ (z >> 31);
}

uint64_t 
crpx_rng_lehmer_seed128 (void *vstate)
{
  __uint128_t state;
  memcpy(&state, vstate, sizeof(state));
  state *= UINT64_C(0xda942042e4dd58b5);
  memcpy(vstate, &state, sizeof(state));
  return (uint64_t)(state >> 64);
}

uint64_t 
crpx_rng_wyrand_seed64 (void *vstate)
{
  uint64_t *state = (uint64_t *) vstate;
  *state += UINT64_C(0xa0761d6478bd642f);
  __uint128_t t = (__uint128_t)(*state) * ((*state) ^ UINT64_C(0xe7037ed1a0b428db));
  return (uint64_t)((t >> 64) ^ t);
}

uint64_t 
crpx_rng_jenkins13_seed256 (void *vstate) // 13 bits of avalanche
{
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
{
  uint64_t *s = (uint64_t *) vstate; 
  uint64_t e = s[0] - ROTL64(s[1], 7);
  s[0] = s[1] ^ ROTL64(s[2], 13);
  s[1] = s[2] + ROTL64(s[3], 37);
  s[2] = s[3] + e;
  s[3] = e + s[0];
  return s[3];
}



uint64_t
crpx_xoroshiro_pv8_seed128 (void *vstate) // 128+ V 2018
{
  uint64_t *v = (uint64_t *) vstate;
  uint64_t v0 = v[0], v1 = v[1];
  uint64_t result = v[0] + v[1]; 
  v1 ^= v0; v[0] = ROTL64(v0, 24) ^ v1 ^ (v1 << 16); v[1] = ROTL64(v1, 37);
  return result;
}

uint64_t
crpx_xoroshiro_pp_seed128 (void *vstate) // 128++
{
  uint64_t *v = (uint64_t *) vstate;
  uint64_t v0 = v[0], v1 = v[1];
  uint64_t result = ROTL64(v0 + v1, 17) + v0; 
  v1 ^= v0; v[0] = ROTL64(v0, 49) ^ v1 ^ (v1 << 21); v[1] = ROTL64(v1, 28);
  return result;
}

uint64_t
crpx_xoroshiro_pp_seed256 (void *vstate) 
{
  uint64_t *v = (uint64_t *) vstate;
  uint64_t t = v[1] << 17; 
  uint64_t result = ROTL64(v[0] + v[3], 23) + v[0];
  v[2] ^= v[0]; v[3] ^= v[1]; v[1] ^= v[2]; v[0] ^= v[3]; v[2] ^= t; v[3] = ROTL64(v[3], 45);
  return result;
}

uint64_t
crpx_xoroshiro_star_seed256 (void *vstate) 
{
  uint64_t *v = (uint64_t *) vstate;
  uint64_t t = v[1] << 17; 
  uint64_t result = ROTL64(v[1] * 5, 7) * 9;
  v[2] ^= v[0]; v[3] ^= v[1]; v[1] ^= v[2]; v[0] ^= v[3]; v[2] ^= t; v[3] = ROTL64(v[3], 45);
  return result;
}

/* Mark Overton, Romu (2020), public-domain source at romu-random.org. State must contain a non-zero word. */
uint64_t
crpx_rng_romu_seed256 (void *vstate)
{ // RomuQuad: rotate the previous fourth word into the first and third updates.
  uint64_t *state = vstate;
  uint64_t x = state[0], y = state[1], z = state[2], w = state[3];
  state[0] = z + ROTL64(w, 52);
  state[1] = y - x;
  state[2] = ROTL64(y + w, 19);
  state[3] = UINT64_C(15241094284759029579) * z;
  return x;
}

uint64_t
crpx_rng_romu_seed192 (void *vstate)
{ // RomuTrio, with three 64-bit state words.
  uint64_t *state = vstate;
  uint64_t x = state[0], y = state[1], z = state[2];
  state[0] = UINT64_C(15241094284759029579) * z;
  state[1] = ROTL64(y - x, 12);
  state[2] = ROTL64(z - y, 44);
  return x;
}

uint64_t
crpx_rng_romu_seed128 (void *vstate)
{ // RomuDuo, with two 64-bit state words.
  uint64_t *state = vstate;
  uint64_t x = state[0], y = state[1];
  state[0] = UINT64_C(15241094284759029579) * y;
  state[1] = ROTL64(y, 36) + ROTL64(y, 15) - x;
  return x;
}

/* Sebastiano Vigna's public-domain xoroshiro128+ (2016) and xorshift recurrences.
   https://prng.di.unimi.it/ ; the 2016 xoroshiro variant is retained for existing explicit selections. */
uint64_t
crpx_xoroshiro_pv6_seed128 (void *vstate)
{
  uint64_t *state = vstate;
  uint64_t first = state[0], second = state[1], result = first + second;
  second ^= first;
  state[0] = ROTL64(first, 55) ^ second ^ (second << 14);
  state[1] = ROTL64(second, 36);
  return result;
}

uint64_t
crpx_rng_xorshift_star_seed64 (void *vstate)
{ // xorshift64*, Vigna (2014); a zero initial state is forbidden.
  uint64_t *state = vstate;
  uint64_t value = *state;
  value ^= value >> 12;
  value ^= value << 25;
  value ^= value >> 27;
  *state = value;
  return value * UINT64_C(2685821657736338717);
}

uint64_t
crpx_rng_xorshift_p_seed128 (void *vstate)
{ // Historical (23,18,5) xorshift128+ variant, with two 64-bit state words.
  uint64_t *state = vstate;
  uint64_t first = state[0], second = state[1];
  first ^= first << 23;
  first ^= first >> 18;
  first ^= second ^ (second >> 5);
  state[0] = second;
  state[1] = first;
  return first + second;
}

#define PCG_DEFAULT_MULTIPLIER_128 ((((__uint128_t)2549297995355413924ULL) << 64) + 4865540595714422341ULL)
uint64_t
crpx_rng_pcg_seed256 (void *vstate)
{ // https://github.com/lemire/testingRNG/blob/master/source/pcg64.h
  __uint128_t s[2];
  memcpy(s, vstate, sizeof(s));
  uint64_t value;
  unsigned int rot;
  s[0] = s[0] * PCG_DEFAULT_MULTIPLIER_128 + s[1]; // s[1] never changes but _must_be_ odd;
  value = ((uint64_t)(s[0] >> 64u)) ^ (uint64_t) s[0];
  rot = s[0] >> 122u;
  memcpy(vstate, s, sizeof(s));
  return (value >> rot) | (value << ((64u - rot) & 63u));
}

uint64_t
crpx_rng_mt19937_seed2504 (void *state) // needs 312 uint64_t for random state and last one is a counter
{ // Nishimura and Matsumoto, MT19937-64 (2004), BSD-3-Clause; see LICENSES/.
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

void 
crpx_rng_pcg_set_seed256 (void *vstate, uint64_t seed)
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
  __uint128_t states[2] = { initstate, initseq };
  memcpy(vstate, states, sizeof(states));
}

void
crpx_rng_mt19937_set_seed2504 (void *state, uint64_t seed)
{ // Nishimura and Matsumoto, MT19937-64 (2004), BSD-3-Clause; see LICENSES/.
  uint64_t *r = (uint64_t *) state;
  if (!seed) seed = 0x8fc18365c966079ULL; // arbitrary prime number by @leomrtns
  r[312] = 313; // counter, forcing generation of 312 words at first call

  for (int i = 0; i < 312; i++) r[i] = crpx_rng_splitmix_seed64 (&seed);
}

void
crpx_rng_xoshiro_jump_seed256 (void *vstate)
{ // Advance by 2^128 draws, using Blackman and Vigna's published jump polynomial.
  static const uint64_t polynomial[4] = {
    UINT64_C(0x180ec6d33cfd0aba), UINT64_C(0xd5a61266f0c9392c),
    UINT64_C(0xa9582618e03fc9aa), UINT64_C(0x39abdc4529b1661c)
  };
  uint64_t *state = vstate, jumped[4] = {0};
  for (size_t word = 0; word < 4; word++) for (unsigned bit = 0; bit < 64; bit++) {
    if (polynomial[word] & (UINT64_C(1) << bit))
      for (size_t i = 0; i < 4; i++) jumped[i] ^= state[i];
    crpx_xoroshiro_pp_seed256 (state);
  }
  memcpy (state, jumped, sizeof (jumped));
}
