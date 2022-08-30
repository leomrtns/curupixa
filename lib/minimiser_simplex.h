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

/*! \file minimiser_simplex.h 
 *  \brief  The Simplex method of Nelder and Mead, also known as the polytope search alogorithm, from the GSL (GPL-3.0)
 * Ref: Nelder, J.A., Mead, R., Computer Journal 7 (1965) pp. 308-313. This implementation uses n+1 corner points in the simplex. 
 * This is the "simplex2" algorithm from GSL, where the size of simplex is calculated as the RMS distance of each vertex from the center 
 * rather than the mean distance. */

#ifndef _curupixa_minimiser_simplex_h_
#define _curupixa_minimiser_simplex_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "global/global_variable.h"

/* Copyright (C) 2007, 2008, 2009 Brian Gough <bjg@network-theory.co.uk>, from the GSL (GPL-3.0)
 * Copyright (C) 2002 Tuomo Keskitalo <tuomo.keskitalo@iki.fi>, Ivo Alxneit  <ivo.alxneit@psi.ch>, from the GSL (GPL-3.0) */

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
