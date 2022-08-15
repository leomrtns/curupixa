/* This test file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include <curupixa.h>
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

// use it like memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
char filename[2048] = "./files/"; // now we can memcpy() file names _after_ prefix_size
const size_t prefix_size = 8; // all modifications to filename[] come after prefix_size


START_TEST(have_instructions)
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

#ifdef HAVE_SSE
    printf("HAVE_SSE from config.h\n");
#endif
#ifdef HAVE_AVX
    printf("HAVE_AVX from config.h\n");
#endif

    printf ("sse   = %d\n", __builtin_cpu_supports("sse"));
    printf ("sse2  = %d\n", __builtin_cpu_supports("sse2"));
    printf ("sse3  = %d\n", __builtin_cpu_supports("sse3"));
    printf ("ssse3 = %d\n", __builtin_cpu_supports("ssse3"));
    printf ("sse4.1= %d\n", __builtin_cpu_supports("sse4.1"));
    printf ("sse4.2= %d\n", __builtin_cpu_supports("sse4.2"));
    printf ("avx   = %d\n", __builtin_cpu_supports("avx"));
    printf ("avx2  = %d\n", __builtin_cpu_supports("avx2"));
  // if something goes wrong you can use: ck_abort_msg ("message");
}
END_TEST

START_TEST (test_mm128)
{
#ifdef __SSE__
    __m128i values =  _mm_setr_epi32(10, 20, 30, 40);
    uint16_t res[8];
    memcpy(res, &values, sizeof(values));
    for (int i = 0; i < 8; i++) {
      printf ("%i mm128 %d\n", i, res[i]);
    }
    printf ("sizeof values is %lu\n", sizeof(values));
    // printf ("extract mm128 0:%d\n", _mm_extract_epi32(values, 0)); // inlining failed with SSE2 but not with AVX2
#endif
}
END_TEST

START_TEST(test_mm256) 
{
#ifdef __AVX2__
    __m256i values2 =  _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8); // set and extract are "int", other functions may need difference between epi and epu (unsigned)
    printf ("0 mm256 %d\n", _mm256_extract_epi32(values2, 0));
    printf ("1 mm256 %d\n", _mm256_extract_epi32(values2, 1));
    printf ("2 mm256 %d\n", _mm256_extract_epi32(values2, 2));
    printf ("3 mm256 %d\n", _mm256_extract_epi32(values2, 3));
    printf ("4 mm256 %d\n", _mm256_extract_epi32(values2, 4));
    printf ("5 mm256 %d\n", _mm256_extract_epi32(values2, 5));
    printf ("6 mm256 %d\n", _mm256_extract_epi32(values2, 6));
    printf ("7 mm256 %d\n", _mm256_extract_epi32(values2, 7));
#endif
}
END_TEST

Suite * this_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("SSE/AVX instructions");
  tc_case = tcase_create("macros");
  tcase_add_test(tc_case, have_instructions);
  tcase_add_test(tc_case, test_mm128);
  tcase_add_test(tc_case, test_mm256);
  suite_add_tcase(s, tc_case);
  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (this_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}
