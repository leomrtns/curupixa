/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file hash_metro.c
    MetroHash 64/128 v1/v2 from J. Andrew Rogers' original Apache-2.0 source (2015-2018).
    Modified for Curupixa: C17, unaligned little-endian reads, native uint64_t seeds and optional 128-bit output.
    Copyright 2015-2018 J. Andrew Rogers. See LICENSES/Apache-2.0.txt and docs/random-hashing.md. */
#include "hash_functions_generators.h"
#include "hash_bytes_internal.h"
#include "internal_random_constants.h"

uint64_t
crpx_metrohash64_v1_seed64 (const void *key, size_t len, const void *seed)
{ // One-shot MetroHash; no allocation or persistent state.
  uint64_t seed_value;
  memcpy (&seed_value, seed, sizeof (seed_value));
  static const uint64_t k0 = 0xC83A91E1;
  static const uint64_t k1 = 0x8648DBDB;
  static const uint64_t k2 = 0x7BDEC03B;
  static const uint64_t k3 = 0x2F5870A5;

  const uint8_t * ptr = (const uint8_t *)(key ? key : "");
  const uint8_t * const end = ptr + len;

  uint64_t hash = ((seed_value + k2) * k0) + len;

  if (len >= 32) {
    uint64_t v[4];
    v[0] = v[1] = v[2] = v[3] = hash;

    do {
      v[0] += crpx_hash_read64 (ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += crpx_hash_read64 (ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += crpx_hash_read64 (ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= end - 32);

    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 33) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 33) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 33) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 33) * k0;
    hash += v[0] ^ v[1];
  }

  if ((end - ptr) >= 16) {
    uint64_t v0 = hash + (crpx_hash_read64 (ptr) * k0); ptr += 8; v0 = ROTR64(v0,33) * k1;
    uint64_t v1 = hash + (crpx_hash_read64 (ptr) * k1); ptr += 8; v1 = ROTR64(v1,33) * k2;
    v0 ^= ROTR64(v0 * k0, 35) + v1;
    v1 ^= ROTR64(v1 * k3, 35) + v0;
    hash += v1;
  }

  if ((end - ptr) >= 8) {
    hash += crpx_hash_read64 (ptr) * k3; ptr += 8;
    hash ^= ROTR64(hash, 33) * k1;

  }

  if ((end - ptr) >= 4) {
    hash += crpx_hash_read32 (ptr) * k3; ptr += 4;
    hash ^= ROTR64(hash, 15) * k1;
  }

  if ((end - ptr) >= 2) {
    hash += crpx_hash_read16 (ptr) * k3; ptr += 2;
    hash ^= ROTR64(hash, 13) * k1;
  }

  if ((end - ptr) >= 1) {
    hash += *ptr * k3;
    hash ^= ROTR64(hash, 25) * k1;
  }

  hash ^= ROTR64(hash, 33);
  hash *= k0;
  hash ^= ROTR64(hash, 33);

  return hash;
}

uint64_t
crpx_metrohash64_v2_seed64 (const void *key, size_t len, const void *seed)
{ // One-shot MetroHash; no allocation or persistent state.
  uint64_t seed_value;
  memcpy (&seed_value, seed, sizeof (seed_value));
  static const uint64_t k0 = 0xD6D018F5;
  static const uint64_t k1 = 0xA2AA033B;
  static const uint64_t k2 = 0x62992FC1;
  static const uint64_t k3 = 0x30BC5B29;

  const uint8_t * ptr = (const uint8_t *)(key ? key : "");
  const uint8_t * const end = ptr + len;

  uint64_t hash = ((seed_value + k2) * k0) + len;

  if (len >= 32) {
    uint64_t v[4];
    v[0] = v[1] = v[2] = v[3] = hash;

    do {
      v[0] += crpx_hash_read64 (ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += crpx_hash_read64 (ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += crpx_hash_read64 (ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= end - 32);

    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 30) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 30) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 30) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 30) * k0;
    hash += v[0] ^ v[1];
  }

  if ((end - ptr) >= 16) {
    uint64_t v0 = hash + (crpx_hash_read64 (ptr) * k2); ptr += 8; v0 = ROTR64(v0,29) * k3;
    uint64_t v1 = hash + (crpx_hash_read64 (ptr) * k2); ptr += 8; v1 = ROTR64(v1,29) * k3;
    v0 ^= ROTR64(v0 * k0, 34) + v1;
    v1 ^= ROTR64(v1 * k3, 34) + v0;
    hash += v1;
  }

  if ((end - ptr) >= 8) {
    hash += crpx_hash_read64 (ptr) * k3; ptr += 8;
    hash ^= ROTR64(hash, 36) * k1;
  }

  if ((end - ptr) >= 4) {
    hash += crpx_hash_read32 (ptr) * k3; ptr += 4;
    hash ^= ROTR64(hash, 15) * k1;
  }

  if ((end - ptr) >= 2) {
    hash += crpx_hash_read16 (ptr) * k3; ptr += 2;
    hash ^= ROTR64(hash, 15) * k1;
  }

  if ((end - ptr) >= 1) {
    hash += *ptr * k3;
    hash ^= ROTR64(hash, 23) * k1;
  }

  hash ^= ROTR64(hash, 28);
  hash *= k0;
  hash ^= ROTR64(hash, 29);

  return hash;
}

uint64_t
crpx_metrohash128_v1_seed64 (const void *key, size_t len, const void *seed, void *out)
{ // One-shot MetroHash; no allocation or persistent state.
  uint64_t seed_value;
  memcpy (&seed_value, seed, sizeof (seed_value));
  static const uint64_t k0 = 0xC83A91E1;
  static const uint64_t k1 = 0x8648DBDB;
  static const uint64_t k2 = 0x7BDEC03B;
  static const uint64_t k3 = 0x2F5870A5;

  const uint8_t * ptr = (const uint8_t *)(key ? key : "");
  const uint8_t * const end = ptr + len;

  uint64_t v[4];

  v[0] = ((seed_value - k0) * k3) + len;
  v[1] = ((seed_value + k1) * k2) + len;

  if (len >= 32) {
    v[2] = ((seed_value + k0) * k2) + len;
    v[3] = ((seed_value - k1) * k3) + len;

    do {
      v[0] += crpx_hash_read64 (ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += crpx_hash_read64 (ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += crpx_hash_read64 (ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= end - 32);

    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 26) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 26) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 26) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 30) * k0;
  }

  if ((end - ptr) >= 16) {
    v[0] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],33) * k3;
    v[1] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[1] = ROTR64(v[1],33) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 17) * k1;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 17) * k0;
  }

  if ((end - ptr) >= 8) {
    v[0] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],33) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 20) * k1;
  }

  if ((end - ptr) >= 4) {
    v[1] += crpx_hash_read32 (ptr) * k2; ptr += 4; v[1] = ROTR64(v[1],33) * k3;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 18) * k0;
  }

  if ((end - ptr) >= 2) {
    v[0] += crpx_hash_read16 (ptr) * k2; ptr += 2; v[0] = ROTR64(v[0],33) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 24) * k1;
  }

  if ((end - ptr) >= 1) {
    v[1] += *ptr * k2; v[1] = ROTR64(v[1],33) * k3;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 24) * k0;
  }

  v[0] += ROTR64((v[0] * k0) + v[1], 13);
  v[1] += ROTR64((v[1] * k1) + v[0], 37);
  v[0] += ROTR64((v[0] * k2) + v[1], 13);
  v[1] += ROTR64((v[1] * k3) + v[0], 37);

  if (out) memcpy (out, v, 16);
  return crpx_mumhash64_mixer (v[0], v[1]);
}

uint64_t
crpx_metrohash128_v2_seed64 (const void *key, size_t len, const void *seed, void *out)
{ // One-shot MetroHash; no allocation or persistent state.
  uint64_t seed_value;
  memcpy (&seed_value, seed, sizeof (seed_value));
  static const uint64_t k0 = 0xD6D018F5;
  static const uint64_t k1 = 0xA2AA033B;
  static const uint64_t k2 = 0x62992FC1;
  static const uint64_t k3 = 0x30BC5B29;

  const uint8_t * ptr = (const uint8_t *)(key ? key : "");
  const uint8_t * const end = ptr + len;

  uint64_t v[4];

  v[0] = ((seed_value - k0) * k3) + len;
  v[1] = ((seed_value + k1) * k2) + len;

  if (len >= 32) {
    v[2] = ((seed_value + k0) * k2) + len;
    v[3] = ((seed_value - k1) * k3) + len;

    do {
      v[0] += crpx_hash_read64 (ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += crpx_hash_read64 (ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += crpx_hash_read64 (ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= end - 32);

    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 33) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 33) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 33) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 33) * k0;
  }

  if ((end - ptr) >= 16) {
    v[0] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],29) * k3;
    v[1] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[1] = ROTR64(v[1],29) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 29) * k1;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 29) * k0;
  }

  if ((end - ptr) >= 8) {
    v[0] += crpx_hash_read64 (ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],29) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 29) * k1;
  }

  if ((end - ptr) >= 4) {
    v[1] += crpx_hash_read32 (ptr) * k2; ptr += 4; v[1] = ROTR64(v[1],29) * k3;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 25) * k0;
  }

  if ((end - ptr) >= 2) {
    v[0] += crpx_hash_read16 (ptr) * k2; ptr += 2; v[0] = ROTR64(v[0],29) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 30) * k1;
  }

  if ((end - ptr) >= 1) {
     v[1] += *ptr * k2; v[1] = ROTR64(v[1],29) * k3;
     v[1] ^= ROTR64((v[1] * k3) + v[0], 18) * k0;
  }

  v[0] += ROTR64((v[0] * k0) + v[1], 33);
  v[1] += ROTR64((v[1] * k1) + v[0], 33);
  v[0] += ROTR64((v[0] * k2) + v[1], 33);
  v[1] += ROTR64((v[1] * k3) + v[0], 33);

  if (out) memcpy (out, v, 16);
  return crpx_mumhash64_mixer (v[0], v[1]);
}
