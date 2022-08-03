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

/*! \file maths_and_bits.h 
 *  \brief Low-level mathematical and bitwise functions.
 *  
 *  This file shold go together with `lowlevel.h` since there are no dependencies. 
 *  Includes a few tricks from https://graphics.stanford.edu/~seander/bithacks.html
 */

#ifndef _curupixa_maths_and_bits_h_
#define _curupixa_maths_and_bits_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "lowlevel.h"

#define CRPX_set_or_clear_bits (w,m,f) (w ^ (-f ^ w) & m) /*!< \brief nonbranching f=set mask bits m: (if (f) w |= m; else w &= ~m)  */
#define CRPX_merge_bits_using_mask (a,b,mask) (a ^ ((a ^ b) & mask)) /*!< \brief mask has 1 where bits from b are selected, and zero where from a: ((a & ~mask) | (b & mask)) */
#define CRPX_reverse_bits_in_byte (b) (((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32) /*!< \brief reverse the 8 bits in a byte (i.e. backwards) using int64 */
#define CRPX_int32_haszerobyte(v) (((v) - 0x01010101UL) & ~(v) & 0x80808080UL) /*!< \brief check if any byte (8bit chunks) in an int32 equals zero */

extern uint32_t crpx_prev_power_of_two (uint32_t x);
extern uint32_t crpx_next_power_of_two (uint32_t x);
extern int crpx_choose_n_k (int n, int k);
void crpx_ordered_combination_n_k (int* result, int n, int k, int order);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
