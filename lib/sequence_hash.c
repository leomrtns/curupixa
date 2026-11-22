/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file sequence_hash.c
    ntHash2 split-rotation recurrence (Kazemi et al., Bioinformatics 38:4812, 2022), generalised to supplied alphabets.
    The one-step rotations and DNA seed constants follow Hamid Mohamadi's MIT-licensed ntHash source; notice retained
    in LICENSES/ntHash-MIT.txt. Changes: C17, explicit iterator state, arbitrary byte maps and seeded protein/text modes. */
#include "sequence_hash.h"
#include "global/hash_functions_generators.h"
#include <stdlib.h>
#include <string.h>

struct crpx_sequence_hash {
  const unsigned char *data;
  size_t window, length, position, valid;
  uint64_t forward, reverse, table[256], outgoing[256], reverse_incoming[4];
  uint16_t codes[256];
  bool dna;
};
static uint64_t split_left (uint64_t value);
static uint64_t split_right (uint64_t value);
static uint64_t split_rotate (uint64_t value, size_t distance);

crpx_status
crpx_sequence_hash_new (size_t window, uint64_t seed, crpx_sequence_alphabet alphabet, crpx_sequence_hash **out)
{ // Tables depend on the iterator, not on hidden global state or the input length.
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!window || alphabet < CRPX_HASH_BYTES || alphabet > CRPX_HASH_PROTEIN) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_sequence_hash *iterator = calloc (1, sizeof (*iterator));
  if (!iterator) return CRPX_ERROR_ALLOCATION;
  iterator->window = window; iterator->dna = alphabet == CRPX_HASH_DNA;
  uint64_t state = seed;
  for (size_t i = 0; i < 256; i++) {
    iterator->table[i] = crpx_hashint_splitmix64 (state += UINT64_C(0x9e3779b97f4a7c15));
    iterator->codes[i] = alphabet == CRPX_HASH_BYTES ? (uint16_t)i : 256;
  }
  if (alphabet == CRPX_HASH_DNA) {
    const unsigned char bases[] = "ACGTU";
    for (size_t i = 0; i < 5; i++) {
      iterator->codes[bases[i]] = iterator->codes[bases[i] + ('a' - 'A')] = (uint16_t)(i == 4 ? 3 : i);
    }
    if (!seed) {
      const uint64_t original[] = {UINT64_C(0x3c8bfbb395c60474), UINT64_C(0x3193c18562a02b4c),
                                   UINT64_C(0x20323ed082572324), UINT64_C(0x295549f54be24456)};
      memcpy (iterator->table, original, sizeof (original));
    }
  } else if (alphabet == CRPX_HASH_PROTEIN) {
    for (unsigned i = 0; i < 26; i++) iterator->codes['A' + i] = iterator->codes['a' + i] = (uint16_t)i;
  }
  for (size_t i = 0; i < 256; i++) iterator->outgoing[i] = split_rotate (iterator->table[i], window);
  if (iterator->dna)
    for (size_t i = 0; i < 4; i++) iterator->reverse_incoming[i] = split_rotate (iterator->table[3 - i], window - 1);
  *out = iterator;
  return CRPX_OK;
}

void
crpx_sequence_hash_free (crpx_sequence_hash *iterator)
{
  free (iterator);
}

crpx_status
crpx_sequence_hash_reset (crpx_sequence_hash *iterator, const void *data, size_t length)
{
  if (!iterator || (!data && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  iterator->data = data; iterator->length = length;
  iterator->position = iterator->valid = 0; iterator->forward = iterator->reverse = 0;
  return CRPX_OK;
}

crpx_status
crpx_sequence_hash_set_map (crpx_sequence_hash *iterator, const uint16_t codes[256])
{ // Equal codes share a hash entry; invalid characters restart the next contiguous window.
  if (!iterator || !codes) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < 256; i++)
    if (codes[i] > 256 || (iterator->dna && codes[i] > 3 && codes[i] != 256)) return CRPX_ERROR_INVALID_ARGUMENT;
  memcpy (iterator->codes, codes, sizeof (iterator->codes));
  return crpx_sequence_hash_reset (iterator, iterator->data, iterator->length);
}

crpx_status
crpx_sequence_hash_next (crpx_sequence_hash *iterator, crpx_sequence_hash_value *value)
{ // Roll both DNA strands together; for protein/text only the forward recurrence is evaluated.
  if (!iterator || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  if (iterator->window > iterator->length) return CRPX_ERROR_NOT_FOUND;
  while (iterator->position < iterator->length) {
    size_t position = iterator->position++;
    unsigned code = iterator->codes[iterator->data[position]];
    if (code == 256) { iterator->valid = 0; iterator->forward = iterator->reverse = 0; continue; }
    if (iterator->valid < iterator->window) {
      iterator->forward = split_left (iterator->forward) ^ iterator->table[code];
      if (iterator->dna) iterator->reverse ^= split_rotate (iterator->table[3 - code], iterator->valid);
      iterator->valid++;
    } else {
      unsigned old = iterator->codes[iterator->data[position - iterator->window]];
      iterator->forward = split_left (iterator->forward) ^ iterator->outgoing[old] ^ iterator->table[code];
      if (iterator->dna)
        iterator->reverse = split_right (iterator->reverse ^ iterator->table[3 - old]) ^ iterator->reverse_incoming[code];
    }
    if (iterator->valid < iterator->window) continue;
    crpx_sequence_hash_value result = {0};
    result.position = iterator->position - iterator->window;
    result.forward = iterator->forward; result.reverse_complement = iterator->reverse;
    result.reverse_is_canonical = iterator->dna && iterator->reverse < iterator->forward;
    result.hash = iterator->dna ? iterator->forward + iterator->reverse : iterator->forward;
    *value = result;
    return CRPX_OK;
  }
  return CRPX_ERROR_NOT_FOUND;
}

static uint64_t
split_left (uint64_t value)
{ // Rotate the lower 33 and upper 31 bits independently, as in ntHash2.
  uint64_t carry = ((value & UINT64_C(0x8000000000000000)) >> 30) | ((value & UINT64_C(0x100000000)) >> 32);
  return ((value << 1) & UINT64_C(0xfffffffdffffffff)) | carry;
}

static uint64_t
split_right (uint64_t value)
{
  uint64_t carry = ((value & UINT64_C(0x200000000)) << 30) | ((value & 1) << 32);
  return ((value >> 1) & UINT64_C(0xfffffffeffffffff)) | carry;
}

static uint64_t
split_rotate (uint64_t value, size_t distance)
{ // Reduce each shift modulo its own subword width; zero shifts never invoke a full-width C shift.
  const uint64_t mask33 = UINT64_C(0x1ffffffff), mask31 = UINT64_C(0x7fffffff);
  uint64_t low = value & mask33, high = value >> 33;
  unsigned left = (unsigned)(distance % 33), right = (unsigned)(distance % 31);
  if (left) low = ((low << left) | (low >> (33 - left))) & mask33;
  if (right) high = ((high << right) | (high >> (31 - right))) & mask31;
  return low | (high << 33);
}
