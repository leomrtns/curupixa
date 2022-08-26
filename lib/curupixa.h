/* 
 *This file is part of curupixa, a low-level library for phylogenomic analysis.
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

/*! \file curupixa.h 
 *  \brief headers exposed to other programs
 */

#ifndef _curupixa_toplevel_h_ 
#define _curupixa_toplevel_h_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "lowlevel.h"
#include "maths_and_bits.h"
#include "hash_functions.h" // includes hash_functions_generators.h
#include "random_number.h"  // includes random_number_generators.h

crpx_global_t crpx_global_init (__attribute__((unused)) uint64_t seed, const char *level_string);
void crpx_global_finalise (crpx_global_t cglob);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
