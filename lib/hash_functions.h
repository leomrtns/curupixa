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

/*! \file hash_functions.h 
 *  \brief headers exposed to other programs
 */

#ifndef _curupixa_hash_functions_h_ 
#define _curupixa_hash_functions_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "random_constants.h"

size_t crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size);
void crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4]);
void crpx_get_time_128bits (uint64_t time[2]);
double crpx_update_elapsed_time_128bits (uint64_t past[2]);

extern uint64_t crpx_fastmix64 (uint64_t x);
extern uint64_t crpx_murmurmix64 (uint64_t k);

uint64_t crpx_hash_pearson_seed2048 (void *vkey, size_t len, const void *vseed); // seed must have >= 256 bytes
uint32_t crpx_hash_pseudocrc32_seed8192 (uint32_t crc, void *vkey, size_t len, const void *vseed); // seed must have >= 1024 bytes (256 x 32bits)
uint32_t crpx_hash_fletcher32 (void *vkey, size_t len); // assumes vkey has pair number of number of bytes (o.w. last byte is lost)
uint32_t crpx_hash_jenkins (void *vkey, size_t len);

uint32_t crpx_hash_jenkins_mailund_seed32 (void *vkey, size_t len, void *vseed);
uint32_t crpx_hash_mailund_seed32 (void *vkey, size_t len, void *vseed);
uint32_t crpx_hash_rotating_seed32 (void *vkey, size_t len, void *vseed);

uint64_t crpx_fasthash64_seed64 (const void *vkey, size_t len, void *vseed);
uint32_t crpx_fnv_hash32 (const void* vkey, size_t len);
uint64_t crpx_fnv_hash64 (const void* vkey, size_t len);
uint32_t crpx_hsieh_hash32_seed32 (const void * vkey, size_t len, void *vseed);

uint64_t crpx_metrohash64_v1_seed64 (const void *vkey, size_t vlen, const void *seed); // 32bits seed cast to 64bits
uint64_t crpx_metrohash64_v2_seed64 (const void *vkey, size_t vlen, const void *seed); // 32bits seed cast to 64bits
void crpx_metrohash128_v1_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout); // 32bits seed cast to 64bits
void crpx_metrohash128_v2_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout); // 32bits seed cast to 64bits
void crpx_siphash_seed128 (const void *in, const size_t inlen, const void *seed, uint8_t *out, const size_t outlen);

/* random numbers (depend on a state which is updated as new numbers are generated) */
uint64_t crpx_rng_wyhash64_state64 (void *vstate);
uint64_t crpx_rng_splitmix64_seed64 (void *vstate);
void cprx_rng_abyssinian_set_seed128 (void *vstate, uint32_t seed);
uint32_t crpx_rng_abyssinian32_seed128 (void *vstate); // 2 x uint64_t 
 
uint64_t crpx_rng_romu_seed256 (void *vstate); // romu_quad: 4 x uint64_t 
uint64_t crpx_rng_romu_seed192 (void *vstate); // romu_trio: 3 x uint64_t 
uint64_t crpx_rng_romu_seed128 (void *vstate); // romu_duo: 2 x uint64_t 
uint64_t crpx_xoro128plus_seed128 (void *vstate);
uint64_t crpx_xs64star_seed64 (void *vstate);
uint64_t crpx_xs128plus_seed128 (void *vstate);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
