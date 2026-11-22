/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
#include "sequence_hash.h"
#include "rolling_hash.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); return 1; } } while (0)

static unsigned char
complement (unsigned char value)
{
  switch (value) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
}

int
main (void)
{
  unsigned char sequence[2300], reversed[2300];
  const size_t windows[] = {1, 2, 3, 7, 31, 32, 33, 64, 65, 1023, 1024};
  for (size_t mode = 0; mode < 3; mode++) for (size_t w = 0; w < sizeof (windows) / sizeof (*windows); w++) {
    size_t window = windows[w];
    for (size_t i = 0; i < sizeof (sequence); i++) {
      if (mode == CRPX_HASH_BYTES) sequence[i] = (unsigned char)(i * 137 + i / 3);
      else if (mode == CRPX_HASH_DNA) sequence[i] = (unsigned char)"ACGTTGCAC"[(i * 7 + i / 3) % 9];
      else sequence[i] = (unsigned char)"ACDEFGHIKLMNPQRSTVWY"[(i * 7 + i / 3) % 20];
    }
    crpx_sequence_hash *rolling = NULL, *reference = NULL;
    CHECK (crpx_sequence_hash_new (window, 987, (crpx_sequence_alphabet)mode, &rolling) == CRPX_OK);
    CHECK (crpx_sequence_hash_new (window, 987, (crpx_sequence_alphabet)mode, &reference) == CRPX_OK);
    CHECK (crpx_sequence_hash_reset (rolling, sequence, sizeof (sequence)) == CRPX_OK);
    for (size_t position = 0; position + window <= sizeof (sequence); position++) {
      crpx_sequence_hash_value actual, expected;
      CHECK (crpx_sequence_hash_next (rolling, &actual) == CRPX_OK && actual.position == position);
      CHECK (crpx_sequence_hash_reset (reference, sequence + position, window) == CRPX_OK);
      CHECK (crpx_sequence_hash_next (reference, &expected) == CRPX_OK);
      CHECK (actual.hash == expected.hash && actual.forward == expected.forward && actual.reverse_complement == expected.reverse_complement);
      if (mode == CRPX_HASH_DNA) {
        for (size_t i = 0; i < window; i++) reversed[i] = complement (sequence[position + window - i - 1]);
        CHECK (crpx_sequence_hash_reset (reference, reversed, window) == CRPX_OK);
        CHECK (crpx_sequence_hash_next (reference, &expected) == CRPX_OK);
        CHECK (actual.hash == expected.hash && actual.forward == expected.reverse_complement);
      }
    }
    crpx_sequence_hash_value unchanged = {.position = 99};
    CHECK (crpx_sequence_hash_next (rolling, &unchanged) == CRPX_ERROR_NOT_FOUND && unchanged.position == 99);
    crpx_sequence_hash_free (rolling); crpx_sequence_hash_free (reference);
  }
  crpx_sequence_hash *iterator = NULL;
  CHECK (crpx_sequence_hash_new (3, 0, CRPX_HASH_DNA, &iterator) == CRPX_OK);
  CHECK (crpx_sequence_hash_reset (iterator, "ACGTNacgu", 9) == CRPX_OK);
  const size_t positions[] = {0, 1, 5, 6};
  uint64_t hashes[4];
  for (size_t i = 0; i < 4; i++) {
    crpx_sequence_hash_value value;
    CHECK (crpx_sequence_hash_next (iterator, &value) == CRPX_OK && value.position == positions[i]);
    hashes[i] = value.hash;
  }
  CHECK (hashes[0] == hashes[1] && hashes[0] == hashes[2] && hashes[0] == hashes[3]);
  crpx_sequence_hash_free (iterator);
  CHECK (crpx_sequence_hash_new (2, 12, CRPX_HASH_PROTEIN, &iterator) == CRPX_OK);
  uint16_t map[256];
  for (size_t i = 0; i < 256; i++) map[i] = 256;
  map['A'] = map['a'] = map['G'] = 0; map['C'] = 1;
  CHECK (crpx_sequence_hash_set_map (iterator, map) == CRPX_OK);
  CHECK (crpx_sequence_hash_reset (iterator, "ACGC", 4) == CRPX_OK);
  crpx_sequence_hash_value first, second;
  CHECK (crpx_sequence_hash_next (iterator, &first) == CRPX_OK);
  CHECK (crpx_sequence_hash_next (iterator, &second) == CRPX_OK);
  CHECK (crpx_sequence_hash_next (iterator, &second) == CRPX_OK && first.hash == second.hash);
  crpx_sequence_hash_free (iterator);
  CHECK (crpx_sequence_hash_new (0, 0, CRPX_HASH_BYTES, &iterator) == CRPX_ERROR_INVALID_ARGUMENT && !iterator);
  return 0;
}
