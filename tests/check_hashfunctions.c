/* This test file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include <curupixa.h>
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif
// use it like memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
char filename[2048] = TEST_FILE_DIR; // now we can memcpy() file names _after_ prefix_size
size_t prefix_size = strlen(TEST_FILE_DIR); // all modifications to filename[] come after prefix_size


START_TEST(small_random_seeds)
{
  size_t i, j;
  uint64_t seed[8];
  uint8_t small[64];
  crpx_global_t cglob = crpx_global_init (0, 0, "debug");

  j = crpx_generate_bytesized_random_seeds (cglob, seed, 64);
  memcpy (small, seed, 64);
  for (i=0; i<j;i++) printf ("%3lu %hhu\n", i, small[i]);
  // if something goes wrong you can use: ck_abort_msg ("message");
  crpx_global_finalise (cglob);
}
END_TEST

START_TEST(big_random_seeds)
{
  size_t i, j;
  uint64_t seed[1000];
  uint8_t small[8000];
  int count = 0;
  crpx_global_t cglob = crpx_global_init (0, 0, "debug");

  j = crpx_generate_bytesized_random_seeds (cglob, seed, 8000);
  memcpy (small, seed, 8000);
  for (i=0; i < j;i++) count += (~small[i]) & 1;
  printf ("number of 8bit zeroes = %d out of 8000\n", count);
  count = 0;
  for (i=0; i < j/8;i++) count += (~seed[i]) & 1;
  printf ("number of 64bit zeroes = %d out of 1000\n", count);
  // if something goes wrong you can use: ck_abort_msg ("message");
  crpx_global_finalise (cglob);
}
END_TEST

Suite * this_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("hash functions");
  tc_case = tcase_create("random seeds");
  tcase_add_test(tc_case, small_random_seeds);
  tcase_add_test(tc_case, big_random_seeds);
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
