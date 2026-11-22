/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file sequence_hash.h
    Split-rotation rolling hashes for DNA, protein and byte windows, complementing rolling_hash.h.
    DNA windows can exceed the exact two-bit k=32 limit of kmerhash.h. */
#ifndef CURUPIXA_SEQUENCE_HASH_H
#define CURUPIXA_SEQUENCE_HASH_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef enum { CRPX_HASH_BYTES = 0, CRPX_HASH_DNA = 1, CRPX_HASH_PROTEIN = 2 } crpx_sequence_alphabet;
typedef struct crpx_sequence_hash crpx_sequence_hash;
typedef struct {
  size_t position;
  uint64_t forward, reverse_complement, hash;
  bool reverse_is_canonical;
} crpx_sequence_hash_value;

/* ntHash2-style 33/31-bit split rotation, O(window) initialisation and O(1) updates, without window-sized storage.
   DNA: ACGT/U, case-insensitive, other bytes break windows; hash is the strand-invariant forward+reverse sum.
   Protein: ASCII letters, case-insensitive, other bytes break windows. Bytes: all 256 values, including zero.
   DNA seed zero uses ntHash 2.4.0's original base constants and canonical hashes. Other seeds and alphabets use
   Curupixa's deterministic tables. The split rotation has period 1023; these are 64-bit non-cryptographic hashes.
   Kazemi et al. (2022), doi:10.1093/bioinformatics/btac564. */
crpx_status crpx_sequence_hash_new (size_t window, uint64_t seed, crpx_sequence_alphabet alphabet, crpx_sequence_hash **out);
void crpx_sequence_hash_free (crpx_sequence_hash *iterator);
crpx_status crpx_sequence_hash_reset (crpx_sequence_hash *iterator, const void *data, size_t length);
/* Override the character map to define reduced alphabets. Codes 0..255 share seed-table entries; 256 breaks a window.
   DNA accepts only 0..3 and 256, with complements 3-code. Resets iteration to the start of the borrowed input. */
crpx_status crpx_sequence_hash_set_map (crpx_sequence_hash *iterator, const uint16_t codes[256]);
/* Borrowed bytes must remain valid and unchanged. NOT_FOUND ends iteration; output is unchanged on failure. */
crpx_status crpx_sequence_hash_next (crpx_sequence_hash *iterator, crpx_sequence_hash_value *value);
#ifdef __cplusplus
}
#endif
#endif
