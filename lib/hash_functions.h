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

size_t crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size);
void crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4]);
void crpx_get_time_128bits (uint64_t time[2]);
double crpx_update_elapsed_time_128bits (uint64_t past[2]);

extern uint64_t crpx_fastmix64 (uint64_t x); /*!< \brief compression, _not_ for RNG */
extern uint64_t crpx_murmurmix64 (uint64_t k); /*!< \brief good dieharder propeties */
extern uint64_t crpx_rrmixer64 (uint64_t x); /*!< \brief good dieharder propeties */
extern uint64_t crpx_moremur64 (uint64_t x); /*!< \brief good dieharder propeties */

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
void crpx_metrohash128_v1_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout); // 32bits seed cast to 64bits
void crpx_metrohash128_v2_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout); // 32bits seed cast to 64bits
void crpx_siphash_seed128 (const void *in, const size_t inlen, const void *seed, uint8_t *out, const size_t outlen); /*!< \brief outlen must be 8 or 16 (64 or 128 bits) */

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */

