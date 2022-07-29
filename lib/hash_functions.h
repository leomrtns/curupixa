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

#include "random_constants.h"

size_t crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size);
void crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4]);
void crpx_get_time_128bits (uint64_t time[2]);
double crpx_update_elapsed_time_128bits (uint64_t past[2]);
uint64_t crpx_wyhash64 (uint64_t *seed); // changes seed state (thus a PRNG) 
uint64_t crpx_splitmix64 (uint64_t *seed); // changes seed state (thus a PRNG) 
uint64_t crpx_fmix64 (uint64_t k);
uint64_t crpx_hash_pearson (void *key, size_t len, const void *seed); // seed must have >= 256 bytes
uint32_t crpx_hash_pseudocrc32 (uint32_t crc, void *key, size_t len, const uint32_t *seed); // seed must have >= 256 elements (of 32bits)
uint32_t crpx_hash_fletcher32 (uint16_t const *data, size_t words);
uint32_t crpx_hash_jenkins (void *key, size_t len);
 
#endif
