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
 *  \brief higher level hash functions (use curupixa_global_t etc.) */

#ifndef _global_hash_functions_h_ 
#define _global_hash_functions_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "hash_functions_generators.h"

size_t crpx_generate_bytesized_random_seeds_from_cpu  (crpx_global_t cglob, void *seed, size_t seed_size);
void   crpx_generate_bytesized_random_seeds_from_seed (crpx_global_t cglob, void *seed, size_t seed_size, uint64_t initial_seed);

void crpx_get_time_128bits (uint64_t time[2]);
double crpx_update_elapsed_time_128bits (uint64_t past[2]);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
