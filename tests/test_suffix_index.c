/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "suffix_tree.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

static int verify(const unsigned char *data, size_t length)
{
  crpx_suffix_index *index;
  CHECK(crpx_suffix_index_new(data, length, &index) == CRPX_OK);
  CHECK(crpx_suffix_index_length(index) == length);
  size_t previous = length;
  for (size_t rank = 0; rank < length; rank++) {
    size_t position;
    CHECK(crpx_suffix_index_position(index, rank, &position) == CRPX_OK);
    CHECK(position < length);
    if (rank) {
      size_t common = length - previous < length - position ? length - previous : length - position;
      int comparison = memcmp(data + previous, data + position, common);
      CHECK(comparison < 0 || (comparison == 0 && length - previous < length - position));
    }
    previous = position;
  }
  for (size_t pattern_length = 1; pattern_length <= 4; pattern_length++) {
    for (unsigned bits = 0; bits < (1u << pattern_length); bits++) {
      unsigned char pattern[4];
      for (size_t i = 0; i < pattern_length; i++) pattern[i] = (bits >> i) & 1;
      size_t expected = 0;
      for (size_t i = 0; pattern_length <= length && i <= length - pattern_length; i++)
        if (!memcmp(data + i, pattern, pattern_length)) expected++;
      size_t first = SIZE_MAX, count = SIZE_MAX;
      crpx_status status = crpx_suffix_index_find(index, pattern, pattern_length, &first, &count);
      if (!expected) CHECK(status == CRPX_ERROR_NOT_FOUND && first == SIZE_MAX && count == SIZE_MAX);
      else {
        CHECK(status == CRPX_OK && count == expected);
        for (size_t i = 0; i < count; i++) {
          size_t position;
          CHECK(crpx_suffix_index_position(index, first + i, &position) == CRPX_OK);
          CHECK(pattern_length <= length - position && !memcmp(data + position, pattern, pattern_length));
        }
      }
    }
  }
  size_t expected_repeat = 0;
  for (size_t a = 0; a <= length; a++) for (size_t b = 0; b <= length; b++) {
    size_t expected = 0, actual;
    while (expected < length - a && expected < length - b && data[a + expected] == data[b + expected]) expected++;
    CHECK(crpx_suffix_index_lcp(index, a, b, &actual) == CRPX_OK && actual == expected);
    if (a != b && expected > expected_repeat) expected_repeat = expected;
  }
  size_t position, repeated;
  CHECK(crpx_suffix_index_longest_repeat(index, &position, &repeated) == CRPX_OK);
  CHECK(repeated == expected_repeat);
  crpx_suffix_index_free(index);
  return 0;
}

int main(void)
{
  unsigned char data[10];
  for (size_t length = 0; length <= 10; length++) for (unsigned bits = 0; bits < (1u << length); bits++) {
    for (size_t i = 0; i < length; i++) data[i] = (bits >> i) & 1;
    CHECK(verify(data, length) == 0);
  }
  const unsigned char text[] = "banana";
  CHECK(verify(text, 6) == 0);
  crpx_suffix_index *index = NULL;
  CHECK(crpx_suffix_index_new(data, SIZE_MAX, &index) == CRPX_ERROR_OVERFLOW && index == NULL);
  return 0;
}
