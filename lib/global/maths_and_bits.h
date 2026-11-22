/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: MIT
 */

/*! \file maths_and_bits.h
    \brief Provide numerical and bit-manipulation primitives shared by other modules.
 */
#ifndef _global_maths_and_bits_h_
#define _global_maths_and_bits_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdint.h>

#define CRPX_set_or_clear_bits (w,m,f) (w ^ (-f ^ w) & m) /*!< \brief nonbranching f=set mask bits m: (if (f) w |= m; else w &= ~m)  */
#define CRPX_merge_bits_using_mask (a,b,mask) (a ^ ((a ^ b) & mask)) /*!< \brief mask has 1 where bits from b are selected, and zero where from a: ((a & ~mask) | (b & mask)) */
#define CRPX_reverse_bits_in_byte (b) (((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32) /*!< \brief reverse the 8 bits in a byte (i.e. backwards) using int64 */
#define CRPX_int32_haszerobyte(v) (((v) - 0x01010101UL) & ~(v) & 0x80808080UL) /*!< \brief check if any byte (8bit chunks) in an int32 equals zero */

extern uint32_t crpx_prev_power_of_two (uint32_t x);
extern uint32_t crpx_next_power_of_two (uint32_t x);
extern uint32_t crpx_reverse_bits32 (uint32_t v);
extern uint64_t crpx_interleave_64bits (uint64_t xylo); /* [y4y3y2y1 x4x3x2x1] --> [y4x4 y3x3 y2x2 y1x1] */
extern uint64_t crpx_deinterleave_64bits (uint64_t interleaved); /* [y4y3y2y1 x4x3x2x1] <-- [y4x4 y3x3 y2x2 y1x1] */
/* Returns zero for invalid inputs or a result above INT_MAX. */
extern int crpx_choose_n_k (int n, int k);
/* Legacy one-based combination elements; invalid arguments leave result unchanged. */
void crpx_ordered_combination_n_k (int* result, int n, int k, int order);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
