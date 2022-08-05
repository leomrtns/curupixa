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
 *  \brief hash functions  */

#ifndef _curupixa_hash_functions_h_ 
#define _curupixa_hash_functions_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "random_constants.h"
#include "maths_and_bits.h"

size_t crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size);
void crpx_get_time_128bits (uint64_t time[2]);
double crpx_update_elapsed_time_128bits (uint64_t past[2]);

extern uint64_t crpx_mumhash64_mixer (uint64_t a, uint64_t b);
extern uint64_t crpx_wyhash64_mixer (uint64_t a, uint64_t b);
extern uint32_t crpx_hash_64_to_32 (uint64_t key);

extern uint64_t crpx_hashint_fastmix64 (uint64_t x); /*!< \brief compression, _not_ for RNG */
extern uint64_t crpx_hashint_murmurmix64 (uint64_t k); /*!< \brief good dieharder propeties */
extern uint64_t crpx_hashint_rrmixer64 (uint64_t x); /*!< \brief good dieharder propeties */
extern uint64_t crpx_hashint_moremur64 (uint64_t x); /*!< \brief good dieharder propeties */
extern uint64_t crpx_hashint_staffordmix13 (uint64_t z);
extern uint64_t crpx_hashint_zixmix64 (uint64_t h);

extern uint32_t crpx_hashint_jenkins32 (uint32_t a);
extern uint32_t crpx_hashint_jenkins32_v2 (uint32_t a);
extern uint32_t crpx_hashint_avalanche (uint32_t a);

uint64_t crpx_hash_pearson_seed2048 (const void *vkey, size_t len, const void *vseed); // seed must have >= 256 bytes
uint32_t crpx_hash_pseudocrc32_seed8192 (const void *vkey, size_t len, const void *vseed, uint32_t crc); // seed >= 1024 bytes (256 x 32bits)
uint32_t crpx_hash_fletcher32 (const void *vkey, size_t len); /*!< \brief _not_for RNG */  // len==pair (o.w. last byte is lost); 
uint32_t crpx_hash_jenkins (const void *vkey, size_t len); /*!< \brief good dieharder propeties */

uint32_t crpx_hash_jenkins_mailund_seed32 (const void *vkey, size_t len, void *vseed); /*!< \brief _not_for RNG */
uint32_t crpx_hash_mailund_seed32 (const void *vkey, size_t len, void *vseed); /*!< \brief good dieharder propeties */
uint32_t crpx_hash_rotating_seed32 (const void *vkey, size_t len, void *vseed); /*!< \brief _not_for RNG */

uint64_t crpx_fasthash64_seed64 (const void *vkey, size_t len, void *vseed); /*!< \brief good dieharder propeties */
uint32_t crpx_fnv_hash32 (const void* vkey, size_t len); /*!< \brief _not_for RNG */
uint64_t crpx_fnv_hash64 (const void* vkey, size_t len); /*!< \brief _not_for RNG */
uint32_t crpx_hsieh_hash32_seed32 (const void * vkey, size_t len, void *vseed); /*!< \brief _not_for RNG (not terrible) */

uint64_t crpx_metrohash64_v1_seed64 (const void *vkey, size_t vlen, const void *seed); /*!< \brief good dieharder propeties */ // 32bits seed cast to 64bits
uint64_t crpx_metrohash64_v2_seed64 (const void *vkey, size_t vlen, const void *seed); /*!< \brief good dieharder propeties */ // 32bits seed cast to 64bits
uint64_t crpx_metrohash128_v1_seed64 (const void *vkey, size_t vlen, const void *seed, void *out); // 32bits seed cast to 64bits; out is 128 bits
uint64_t crpx_metrohash128_v2_seed64 (const void *vkey, size_t vlen, const void *seed, void *out); // 32bits seed cast to 64bits; out is 128 bits
uint64_t crpx_murmurhash3_128bits (const void *key, const size_t len, const uint32_t seed, void *out); // out[] is 128 bits
uint32_t crpx_murmurhash3_32bits (const void *data, const size_t nbytes, const uint32_t seed);
uint64_t crpx_siphash128_seed128 (const void *in, const size_t inlen, const void *seed, void *out);
uint64_t crpx_siphash64_seed128 (const void *in, const size_t inlen, const void *seed);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
