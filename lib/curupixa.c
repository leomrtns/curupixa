/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.c 
 *  \brief Lowest level basic functions, that should be available to all other modules. 
 */

#include "curupixa.h"

#ifdef __SSE__
void test_mm128(void);
#endif
#ifdef __AVX2__
void test_mm256(void);
#endif

int main()
{
#ifdef __SSE__
    printf("SSE defined\n");
#endif
#ifdef __SSE2__
    printf("SSE2 defined\n");
#endif
#ifdef __AVX__
    printf("AVX defined\n");
#endif
#ifdef __AVX2__
    printf("AVX2 defined\n");
#endif

#ifdef HAVE_SSE2
    printf("HAVE_SSE2\n");
#endif
#ifdef HAVE_AVX2
    printf("HAVE_AVX2\n");
#endif

    printf ("sse2  = %d\n", __builtin_cpu_supports("sse2"));
    printf ("avx2  = %d\n", __builtin_cpu_supports("avx2"));
#ifdef __SSE__
    if (__builtin_cpu_supports("sse")) test_mm128(); // compiled with SSE4 _and_ runtime host supports it 
#endif
#ifdef __AVX2__
    if (__builtin_cpu_supports("avx2")) test_mm256();
#endif
return 0;
}

#ifdef __SSE__
void test_mm128(void) {
    __m128i values =  _mm_setr_epi32(10, 20, 30, 40);
    uint16_t res[4];
    memcpy(res, &values, sizeof(values));
    for (int i = 0; i < 8; i++) {
      printf ("%i mm128 %d\n", i, res[i]);
    }
    // printf ("extract mm128 0:%d\n", _mm_extract_epi32(values, 0)); // inlining failed with SSE2 but not with AVX2
    return;
}
#endif

#ifdef __AVX2__
void test_mm256(void) {
    __m256i values2 =  _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8); // set and extract are "int", other functions may need difference between epi and epu (unsigned)
    int res;
    for (int i = 0; i < 8; i++) {
      res = _mm256_extract_epi32(values2, (const int) i);
      printf ("%i mm256 %d\n", i, res);
    }
    printf ("extract mm256 0:%d\n", _mm256_extract_epi32(values2, 0));
    return;
}
#endif
