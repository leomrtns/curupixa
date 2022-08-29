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

/*! \file quasi_random.h 
 *  \brief Quasi-random number generator using Halton sequences.
 *  algorithms based on the GSL library (GPL3) */ 

#ifndef _curupixa_quasi_random_h_
#define _curupixa_quasi_random_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "global/global_variable.h"

typedef struct {
  size_t size, leap, rem, iteration; /*!< leap and remainder are used by halton */
  double *r, *ko; /*!< r = vector with quasi-random, ko = korobov aux vector with weights (random generators) */
  crpx_global_t cglob;
} crpx_quasi_random_struct, *crpx_quasi_random_t;

crpx_quasi_random_t new_crpx_quasi_random (crpx_global_t cglob, size_t size);
void del_crpx_quasi_random (crpx_quasi_random_t q);
void crpx_quasi_random_reset (crpx_quasi_random_t q);
void crpx_quasi_random_next_korobov (crpx_quasi_random_t q);
void crpx_quasi_random_next_halton (crpx_quasi_random_t q);
void crpx_quasi_random_next_halton_original (crpx_quasi_random_t q);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
