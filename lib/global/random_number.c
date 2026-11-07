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

#include "global_variable.h"
#include "random_number.h"
#include "runtime_internal.h"

crpx_status
crpx_set_random_generator_checked (crpx_global_t cglob, uint8_t rng_id, uint64_t seed)
{
  if (!cglob || rng_id > 20) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_global_struct configuration = {0};
  switch (rng_id) {
    case 0:  configuration.rng_get = &crpx_rng_wyhash_state64;       configuration.rng_size = 1; strcpy (configuration.rng_name, "0.wyhash_64"); break;
    case 1:  configuration.rng_get = &crpx_rng_lehmer_seed128;       configuration.rng_size = 2; strcpy (configuration.rng_name, "1.lehmer_64"); break;
    case 2:  configuration.rng_get = &crpx_rng_moremur_seed64;       configuration.rng_size = 1; strcpy (configuration.rng_name, "2.moremur_64"); break;
    case 3:  configuration.rng_get = &crpx_rng_splitmix_seed64;      configuration.rng_size = 1; strcpy (configuration.rng_name, "3.splitmix_64"); break;
    case 4:  configuration.rng_get = &crpx_rng_romu_seed128;         configuration.rng_size = 2; strcpy (configuration.rng_name, "4.romu_128"); break;
    case 5:  configuration.rng_get = &crpx_rng_jenkins13_seed256;    configuration.rng_size = 4; strcpy (configuration.rng_name, "5.jenkins13_256"); break;
    case 6:  configuration.rng_get = &crpx_rng_jenkins19_seed256;    configuration.rng_size = 4; strcpy (configuration.rng_name, "6.jenkins19_256"); break;
    case 7:  configuration.rng_get = &crpx_rng_xorshift_star_seed64; configuration.rng_size = 1; strcpy (configuration.rng_name, "7.xorshift_s_64"); break;
    case 8:  configuration.rng_get = &crpx_rng_romu_seed192;         configuration.rng_size = 3; strcpy (configuration.rng_name, "8.romu_192"); break;
    case 9:  configuration.rng_get = &crpx_xoroshiro_pv6_seed128;    configuration.rng_size = 2; strcpy (configuration.rng_name, "9.xoroshiro_pv6_128"); break;
    case 10: configuration.rng_get = &crpx_xoroshiro_pv8_seed128;    configuration.rng_size = 2; strcpy (configuration.rng_name, "10.xoroshiro_pv8_128"); break;
    case 11: configuration.rng_get = &crpx_rng_romu_seed256;         configuration.rng_size = 4; strcpy (configuration.rng_name, "11.romu_256"); break;
    case 12: configuration.rng_get = &crpx_rng_xorshift_p_seed128;   configuration.rng_size = 2; strcpy (configuration.rng_name, "12.xorshift_p_128"); break;
    case 13: configuration.rng_get = &crpx_rng_rrmixer_seed64;       configuration.rng_size = 1; strcpy (configuration.rng_name, "13.rrmixer_64"); break;
    case 14: configuration.rng_get = &crpx_xoroshiro_pp_seed128;     configuration.rng_size = 2; strcpy (configuration.rng_name, "14.xoroshiro_pp_128"); break;
    case 15: configuration.rng_get = &crpx_xoroshiro_star_seed256;   configuration.rng_size = 4; strcpy (configuration.rng_name, "15.xoroshiro_s_256"); break;
    case 16: configuration.rng_get = &crpx_rng_wyrand_seed64;        configuration.rng_size = 1; strcpy (configuration.rng_name, "16.wyrand_64"); break;
    case 17: configuration.rng_get = &crpx_xoroshiro_pp_seed256;     configuration.rng_size = 4; strcpy (configuration.rng_name, "17.xoroshiro_pp_256"); break;
    case 18: configuration.rng_get = &crpx_rng_pcg_seed256;          configuration.rng_size = 4; strcpy (configuration.rng_name, "18.pcg_256"); break;
    case 19: configuration.rng_get = &crpx_rng_xorshift_seed528;    configuration.rng_size = 66; strcpy (configuration.rng_name, "19.xorshift_528"); break;
    default: configuration.rng_get = &crpx_rng_mt19937_seed2504;   configuration.rng_size = 313; strcpy (configuration.rng_name, "20.mt19937"); break;
  }
  size_t bytes = configuration.rng_size * sizeof(uint64_t);
  uint64_t *state = malloc(bytes);
  if (!state) return CRPX_ERROR_ALLOCATION;
  if (!seed && crpx_generate_bytesized_random_seeds_from_cpu(cglob, &seed, sizeof(seed)) != sizeof(seed))
    seed = UINT64_C(0x1576359c7fcbd9df);
  if (rng_id == 20) crpx_rng_mt19937_set_seed2504(state, seed);
  else if (rng_id == 19) crpx_rng_xorshift_set_seed528(state, seed);
  else if (rng_id == 18) crpx_rng_pcg_set_seed256(state, seed);
  else {
    for (size_t i = 0; i < configuration.rng_size; i++) state[i] = crpx_rng_splitmix_seed64(&seed) | 1;
  }
  for (size_t i = 0; i < configuration.rng_size; i++) configuration.rng_get(state);
  crpx_runtime_lock(&cglob->runtime_state->rng_lock);
  uint64_t *previous = cglob->rng_seed_vector;
  cglob->rng_seed_vector = state;
  cglob->rng_get = configuration.rng_get;
  cglob->rng_size = configuration.rng_size;
  memcpy(cglob->rng_name, configuration.rng_name, sizeof(cglob->rng_name));
  free(previous);
  crpx_runtime_unlock(&cglob->runtime_state->rng_lock);
  return CRPX_OK;
}

void
crpx_set_random_generator (crpx_global_t context, uint8_t generator, uint64_t seed)
{
  crpx_status status = crpx_set_random_generator_checked(context, generator, seed);
  if (status != CRPX_OK) crpx_logger_error(context, "Could not set random generator (%d)", (int)status);
}

uint64_t
crpx_random_64bits (crpx_global_t cglob)
{
  if (!cglob) return 0;
  crpx_runtime_lock(&cglob->runtime_state->rng_lock);
  uint64_t value = cglob->rng_get && cglob->rng_seed_vector ? cglob->rng_get(cglob->rng_seed_vector) : 0;
  crpx_runtime_unlock(&cglob->runtime_state->rng_lock);
  return value;
}

uint32_t
crpx_random_32bits (crpx_global_t cglob)
{
  uint64_t h = crpx_random_64bits (cglob);
  return (uint32_t)(0xffffffff & (h - (h >> 32))); // Fermat residue (https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c)
}

uint32_t
crpx_random_32bits_extra (crpx_global_t cglob, uint32_t *extra_result)
{
  uint64_t h = crpx_random_64bits (cglob);
  if (extra_result) *extra_result = h >> 32;
  return (uint32_t)(0xffffffff & h);
}

uint64_t
crpx_random_range (crpx_global_t cglob, uint64_t n)
{
  uint64_t result = 0;
  if (crpx_random_range_checked(cglob, n, &result) != CRPX_OK)
    crpx_logger_error(cglob, "Random range requires a positive bound and a live context");
  return result;
}

/* in https://lemire.me/blog/2017/02/28/how-many-floating-point-numbers-are-in-the-interval-01/ he mentions 
 *    mantissa = 52 (23 for 32bits)  -->  random(0,2^53)/2^53 = [0,1) --> rnd() >> 11 / 9007199254740992
 *    and that rnd(0,1) * N is uniform(0,N) only for N < 2^53   */   
double
crpx_random_double (crpx_global_t cglob) // [0,1)
{
  return (double)(crpx_random_64bits(cglob) >> 11) / 9007199254740992.0; // 2^53
}

double
crpx_random_double_include_one (crpx_global_t cglob) // [0,1]
{
  return (double)(crpx_random_64bits(cglob) >> 11) / 9007199254740991.0; // 2^53-1
}

double
crpx_random_double_positive (crpx_global_t cglob) // (0,1)
{
  if (!cglob) return NAN;
  double x;
  do { x = (crpx_random_64bits(cglob) >> 11) / 9007199254740992.0; } while (x < 2 * DBL_MIN);
 // return ((double)(genrand64_int64(context) >> 12) + 0.5) / 4503599627370496.0; // alternative 
 return x;
}

double
crpx_random_double_positive_include_one (crpx_global_t cglob) // (0,1]
{
  if (!cglob) return NAN;
  double x;
  do { x = (crpx_random_64bits(cglob) >> 11) / 9007199254740991.0; } while (x < 2 * DBL_MIN);
 return x;
}

double
crpx_random_normal (crpx_global_t cglob, double *extra_result) // generates _2_ random numbers, one is stored in extra_result 
{ // from biomcmc
  if (!cglob) return NAN;
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    /*    u and v are between -1 and +1 => -1 + 2 * U       */
    u = -1. + 2. * ((double)(crpx_random_64bits (cglob) >> 11) / 9007199254740991.0); /* (2^53) - 1 => 1 included */
    v = -1. + 2. * ((double)(crpx_random_64bits (cglob) >> 11) / 9007199254740991.0); /* (2^53) - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  if (extra_result) *extra_result = u * s;
  return v * s;
}

double
crpx_random_normal_fast (crpx_global_t cglob, double *extra_result) // generates _2_ random numbers, one is stored in extra_result 
{ // from biomcmc (using 32 bits instead of 53)
  if (!cglob) return NAN;
  double s, u, v; /* u and v are U(-1,1) = -1 + 2 * U(0,1) */
  uint32_t res2;
  do { /* Marsaglia's Polar method; runs, on average, 1.2732 times */
    /*    u and v are between -1 and +1 => -1 + 2 * U       */
    u = -1. + 2. * ((double)(crpx_random_32bits_extra (cglob, &res2)) / 4294967295.0); /* (2^32) - 1 => 1 included */
    v = -1. + 2. * ((double)(res2) / 4294967295.0); /* (2^32) - 1 => 1 included */
    s = (u*u) + (v*v);
  } while ((s <= 0.) || (s >= 1.));
  s = sqrt (-2. * log (s)/s);
  if (extra_result) *extra_result = u * s;
  return v * s;
}


crpx_status
crpx_random_range_checked (crpx_global_t context, uint64_t bound, uint64_t *result)
{
  if (!context || !bound || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t threshold = (uint64_t)(-bound) % bound;
  uint64_t value;
  do { value = crpx_random_64bits(context); } while (value < threshold);
  *result = value % bound;
  return CRPX_OK;
}

void
crpx_rng_stream_init (crpx_rng_stream *stream, uint64_t seed, uint64_t stream_id)
{
  if (!stream) return;
  /* Hash the worker/rank identifier before combining it with the user seed. */
  uint64_t identifier = stream_id;
  stream->state = seed ^ crpx_rng_splitmix_seed64(&identifier);
}

uint64_t
crpx_rng_stream_next (crpx_rng_stream *stream)
{
  return stream ? crpx_rng_splitmix_seed64(&stream->state) : 0;
}

crpx_status
crpx_rng_stream_range (crpx_rng_stream *stream, uint64_t bound, uint64_t *result)
{
  if (!stream || !bound || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t threshold = (uint64_t)(-bound) % bound;
  uint64_t value;
  do { value = crpx_rng_stream_next(stream); } while (value < threshold);
  *result = value % bound;
  return CRPX_OK;
}

double
crpx_rng_stream_double (crpx_rng_stream *stream)
{
  return (double)(crpx_rng_stream_next(stream) >> 11) * 0x1.0p-53;
}
