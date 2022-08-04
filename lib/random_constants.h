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

/*! \file random_constants.h 
 *  \brief vectors of prime numbers and random numbers used by hashes and PRNGs. Examples include the rolling
 *  hash (DNA bases mapped to a random value) and the spice (initial states used to generate several streams)
 *  These functions and macros are _not_ exposed to the user. Every calling function must include this file. */

#ifndef _curupixa_random_constants_h_
#define _curupixa_random_constants_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "lowlevel.h"

#define ROTL64(x, b) (uint64_t)(((x) << (b)) | ((x) >> (64 - (b))))
#define ROTR64(x, b) (uint64_t)(((x) >> (b)) | ((x) << (64 - (b))))
#define ROTL32(x, b) (uint32_t)(((x) << (b)) | ((x) >> (32 - (b))))
#define MIX32(a,b,c) { a -= b; a -= c; a ^= (c>>13);        \
b -= c; b -= a; b ^= (a<<8);  c -= a; c -= b; c ^= (b>>13); \
a -= b; a -= c; a ^= (c>>12); b -= c; b -= a; b ^= (a<<16); \
c -= a; c -= b; c ^= (b>>5);  a -= b; a -= c; a ^= (c>>3);  \
b -= c; b -= a; b ^= (a<<10); c -= a; c -= b; c ^= (b>>15); }

extern uint16_t crpx_random_prime32_length;
extern uint32_t crpx_random_prime32[];
extern uint16_t crpx_random_prime64_length;
extern uint64_t crpx_random_prime64[];
extern uint16_t crpx_random64_length;
extern uint64_t crpx_random64[];

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
