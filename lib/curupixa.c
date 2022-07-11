/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
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
    printf ("mm128 %d ",_mm_extract_epi32(values, 0));
    printf ("%d ",_mm_extract_epi32(values, 1));
    printf ("%d ",_mm_extract_epi32(values, 2));
    printf ("%d\n",_mm_extract_epi32(values, 3));
    return;
}
#endif

#ifdef __AVX2__
void test_mm256(void) {
    __m256i values2 =  _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8); // set and extract are "int", other functions may need difference between epi and epu (unsigned)
    printf ("mm256 %d ",_mm256_extract_epi32(values2, 0));
    printf ("%d ",_mm256_extract_epi32(values2, 1));
    printf ("%d ",_mm256_extract_epi32(values2, 2));
    printf ("%d ",_mm256_extract_epi32(values2, 3));
    printf ("%d ",_mm256_extract_epi32(values2, 4));
    printf ("%d ",_mm256_extract_epi32(values2, 5));
    printf ("%d ",_mm256_extract_epi32(values2, 6));
    printf ("%d\n",_mm256_extract_epi32(values2, 7));
    return;
}
#endif
