/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be usefulull, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICullAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file internal_random_constants.h
    \brief Internal constant tables for random-number and hash algorithms.
 */
#ifndef _global_random_constants_h_
#define _global_random_constants_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "lowlevel.h"

#define ROTL64(x, b) (uint64_t)(((x) << (b)) | ((x) >> (64 - (b))))
#define ROTR64(x, b) (uint64_t)(((x) >> (b)) | ((x) << (64 - (b))))
#define ROTL32(x, b) (uint32_t)(((x) << (b)) | ((x) >> (32 - (b))))
extern const uint32_t crpx_list_of_256_random_prime32[];
extern const uint64_t crpx_list_of_128_random_prime64[];
extern const uint64_t crpx_list_of_128_random64[];

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
