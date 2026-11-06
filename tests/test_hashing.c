/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "hll.h"
#include "kmerhash.h"
#include "rolling_hash.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

static int test_hll(void)
{
  unsigned char sentinel;
  crpx_hll *invalid = (crpx_hll *)&sentinel;
  CHECK(crpx_hll_new(3, 0, &invalid) == CRPX_ERROR_INVALID_ARGUMENT && !invalid);
  invalid = (crpx_hll *)&sentinel;
  CHECK(crpx_hll_deserialize("bad", 3, &invalid) == CRPX_ERROR_PARSE && !invalid);
  crpx_hll *all, *even, *odd;
  CHECK(crpx_hll_new(14, 10, &all) == CRPX_OK);
  CHECK(crpx_hll_new(14, 10, &even) == CRPX_OK);
  CHECK(crpx_hll_new(14, 10, &odd) == CRPX_OK);
  double estimate = -1, error = -1;
  CHECK(crpx_hll_estimate(all, &estimate, &error) == CRPX_OK && estimate == 0);
  CHECK(error > 0 && error < 0.01);
  for (unsigned i = 0; i < 100000; i++) {
    char key[32];
    int length = snprintf(key, sizeof(key), "%u", i);
    CHECK(crpx_hll_add(all, key, (size_t)length) == CRPX_OK);
    CHECK(crpx_hll_add(i & 1 ? odd : even, key, (size_t)length) == CRPX_OK);
  }
  CHECK(crpx_hll_estimate(all, &estimate, NULL) == CRPX_OK);
  CHECK(fabs(estimate - 100000) < 5000);
  CHECK(crpx_hll_merge(even, odd) == CRPX_OK);
  size_t length = 0, merged_length = 0;
  CHECK(crpx_hll_serialize(all, NULL, 0, &length) == CRPX_ERROR_OUT_OF_RANGE);
  unsigned char *reference = malloc(length), *merged = malloc(length);
  CHECK(reference && merged);
  CHECK(crpx_hll_serialize(all, reference, length, &length) == CRPX_OK);
  CHECK(crpx_hll_serialize(even, merged, length, &merged_length) == CRPX_OK);
  CHECK(length == merged_length && !memcmp(reference, merged, length));
  crpx_hll *copy;
  CHECK(crpx_hll_deserialize(reference, length, &copy) == CRPX_OK);
  CHECK(crpx_hll_serialize(copy, merged, length, &merged_length) == CRPX_OK);
  CHECK(!memcmp(reference, merged, length));
  CHECK(crpx_hll_add(copy, "99", 2) == CRPX_OK);
  CHECK(crpx_hll_serialize(copy, merged, length, &merged_length) == CRPX_OK);
  CHECK(!memcmp(reference, merged, length));
  crpx_hll_free(copy);
  copy = NULL;
  reference[24] = 255;
  CHECK(crpx_hll_deserialize(reference, length, &copy) == CRPX_ERROR_PARSE && copy == NULL);
  CHECK(crpx_hll_new(14, 11, &copy) == CRPX_OK);
  CHECK(crpx_hll_merge(even, copy) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_hll_free(copy);
  crpx_hll_reset(all);
  CHECK(crpx_hll_estimate(all, &estimate, NULL) == CRPX_OK && estimate == 0);
  crpx_hll_free(all); crpx_hll_free(even); crpx_hll_free(odd);
  free(reference); free(merged);
  return 0;
}

static unsigned code(unsigned char base)
{
  switch (base) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': case 'U': case 'u': return 3;
    default: return 4;
  }
}

static int test_kmers(void)
{
  const char sequence[] = "acGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNACGTUACG-ACGTACGTACGTACGTACGTACGTACGTACGT";
  const size_t length = sizeof(sequence) - 1;
  for (size_t k = 1; k <= 32; k++) {
    crpx_kmer_iterator *iterator;
    CHECK(crpx_kmer_iterator_new(k, 0, &iterator) == CRPX_OK);
    CHECK(crpx_kmer_iterator_reset(iterator, sequence, length) == CRPX_OK);
    crpx_kmer actual;
    for (size_t position = 0; position + k <= length; position++) {
      uint64_t forward = 0, reverse = 0;
      int valid = 1;
      for (size_t j = 0; j < k; j++) {
        unsigned base = code(sequence[position + j]);
        if (base > 3) { valid = 0; break; }
        forward = (forward << 2) | base;
        reverse |= (uint64_t)(3 - base) << (2 * j);
      }
      if (!valid) continue;
      CHECK(crpx_kmer_iterator_next(iterator, &actual) == CRPX_OK);
      CHECK(actual.position == position && actual.forward == forward && actual.reverse_complement == reverse);
      CHECK(actual.canonical == (forward < reverse ? forward : reverse));
      char reverse_sequence[33];
      const char *alphabet = "ACGT";
      for (size_t j = 0; j < k; j++) reverse_sequence[j] = alphabet[3 - code(sequence[position + k - j - 1])];
      crpx_kmer_iterator *reverse_iterator;
      CHECK(crpx_kmer_iterator_new(k, 0, &reverse_iterator) == CRPX_OK);
      CHECK(crpx_kmer_iterator_reset(reverse_iterator, reverse_sequence, k) == CRPX_OK);
      crpx_kmer counterpart;
      CHECK(crpx_kmer_iterator_next(reverse_iterator, &counterpart) == CRPX_OK);
      CHECK(counterpart.hash == actual.hash && counterpart.canonical == actual.canonical);
      crpx_kmer_iterator_free(reverse_iterator);
    }
    CHECK(crpx_kmer_iterator_next(iterator, &actual) == CRPX_ERROR_NOT_FOUND);
    CHECK(crpx_kmer_iterator_reset(iterator, NULL, 0) == CRPX_OK);
    CHECK(crpx_kmer_iterator_next(iterator, &actual) == CRPX_ERROR_NOT_FOUND);
    crpx_kmer_iterator_free(iterator);
  }
  return 0;
}

static int test_rolling(void)
{
  unsigned char input[300];
  for (size_t i = 0; i < sizeof(input); i++) input[i] = (unsigned char)(i * 13);
  const size_t windows[] = { 1, 2, 31, 32, 63, 64, 65, 128, 300, 301 };
  for (size_t test = 0; test < sizeof(windows) / sizeof(windows[0]); test++) {
    size_t window = windows[test];
    crpx_rolling_hash *iterator, *reference;
    CHECK(crpx_rolling_hash_new(window, 77, &iterator) == CRPX_OK);
    CHECK(crpx_rolling_hash_new(window, 77, &reference) == CRPX_OK);
    CHECK(crpx_rolling_hash_reset(iterator, input, sizeof(input)) == CRPX_OK);
    size_t actual_position, reference_position;
    uint64_t actual, expected;
    for (size_t position = 0; window <= sizeof(input) && position <= sizeof(input) - window; position++) {
      CHECK(crpx_rolling_hash_next(iterator, &actual, &actual_position) == CRPX_OK);
      CHECK(actual_position == position);
      CHECK(crpx_rolling_hash_reset(reference, input + position, window) == CRPX_OK);
      CHECK(crpx_rolling_hash_next(reference, &expected, &reference_position) == CRPX_OK);
      CHECK(actual == expected);
    }
    CHECK(crpx_rolling_hash_next(iterator, &actual, &actual_position) == CRPX_ERROR_NOT_FOUND);
    crpx_rolling_hash_free(iterator);
    crpx_rolling_hash_free(reference);
  }
  return 0;
}

int main(void)
{
  CHECK(test_hll() == 0);
  CHECK(test_kmers() == 0);
  CHECK(test_rolling() == 0);
  return 0;
}
