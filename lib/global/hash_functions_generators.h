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

/*! \file hash_functions_generators.h 
 *  \brief lower level, simple hash functions  */

#ifndef _global_hash_functions_generators_h_ 
#define _global_hash_functions_generators_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "maths_and_bits.h"

extern uint64_t crpx_mumhash64_mixer (uint64_t a, uint64_t b);
extern uint64_t crpx_wyhash64_mixer (uint64_t a, uint64_t b);
extern uint32_t crpx_hash_64_to_32 (uint64_t key);

extern uint64_t crpx_hashint_staffordmix64 (uint64_t z); // same as hashint_splitmix64 but adds prime number as initial state
extern uint64_t crpx_hashint_splitmix64 (uint64_t x); // same as rng_splitmix with state=0 and hashint_staffordmix without state
extern uint64_t crpx_hashint_spitmix64_inverse (uint64_t x);
extern uint64_t crpx_hashint_degski64 (uint64_t x);
extern uint64_t crpx_hashint_degski64_inverse (uint64_t x);
extern uint64_t crpx_hashint_fastmix64 (uint64_t x); /*!< \brief compression, _not_ for RNG */
extern uint64_t crpx_hashint_murmurmix64 (uint64_t k); 
extern uint64_t crpx_hashint_rrmixer64 (uint64_t x); /*!< \brief good dieharder/practrand propeties */
extern uint64_t crpx_hashint_nasam64 (uint64_t x);   /*!< \brief good dieharder/practrand propeties */
extern uint64_t crpx_hashint_pelican64 (uint64_t z); /*!< \brief good dieharder/practrand propeties */
extern uint64_t crpx_hashint_moremur64 (uint64_t x); /*!< \brief good dieharder/practrand propeties */
extern uint64_t crpx_hashint_entropy (uint64_t x);

extern uint32_t crpx_hashint_jenkins (uint32_t a);
extern uint32_t crpx_hashint_jenkins_v2 (uint32_t a);
extern uint32_t crpx_hashint_avalanche (uint32_t a);
extern uint32_t crpx_hashint_murmurmix (uint32_t x);
extern uint32_t crpx_hashint_wellons3ple (uint32_t x); /*!< \brief OK dieharder/practrand propeties (not great) */
extern uint32_t crpx_hashint_wellons3ple_inverse (uint32_t x); // inverse of crpx_hashint_wellons3ple()
extern uint32_t crpx_hashint_wellons (uint32_t x);     /*!< \brief OK dieharder/practrand propeties (not great) */
extern uint32_t crpx_hashint_wellons_inverse (uint32_t x);
extern uint32_t crpx_hashint_degski (uint32_t x);
extern uint32_t crpx_hashint_degski_inverse (uint32_t x);

extern uint16_t crpx_hashint_2xor_16bits (uint16_t x); // 2-round xorshift-multiply; bias = 0.00859
extern uint16_t crpx_hashint_3xor_16bits (uint16_t x); // 3-round xorshift-multiply; bias = 0.00459
extern uint16_t crpx_hashint_noxor_16bits (uint16_t x); // No multiplication; bias = 0.02384

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
uint64_t crpx_siphash128_seed128 (const void *in, const size_t inlen, const void *seed, void *out); // return 64 bits is a mixer of the 128bits, for true 64 bits use siphash64
uint64_t crpx_siphash64_seed128 (const void *in, const size_t inlen, const void *seed);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
