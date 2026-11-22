/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT AND BSD-3-Clause AND CC0-1.0 */
/*! \file hash_bytes.c
    Byte-array hashes with explicit seeds. See docs/random-hashing.md for original sources, licences and changed outputs.
    Tan fast-hash: MIT; Hsieh SuperFastHash: BSD-3-Clause; Jenkins, FNV, MurmurHash and SipHash: public domain/CC0. */
#include "hash_functions_generators.h"
#include "hash_bytes_internal.h"
#include "internal_random_constants.h"

static uint32_t jenkins_bytes (const uint8_t *key, size_t length, uint32_t hash);
static uint32_t jenkins_lookup2 (const uint8_t *key, size_t length, uint32_t hash);

uint64_t
crpx_fnv_hash64 (const void *key, size_t length)
{ // FNV-1: multiply then XOR. The original FNV reference implementation is CC0.
  const uint8_t *bytes = key;
  uint64_t hash = UINT64_C(14695981039346656037);
  for (size_t i = 0; i < length; i++) hash = hash * UINT64_C(1099511628211) ^ bytes[i];
  return hash;
}

uint32_t
crpx_fnv_hash32 (const void *key, size_t length)
{
  const uint8_t *bytes = key;
  uint32_t hash = UINT32_C(2166136261);
  for (size_t i = 0; i < length; i++) hash = hash * UINT32_C(16777619) ^ bytes[i];
  return hash;
}

uint64_t
crpx_fasthash64_seed64 (const void *key, size_t length, void *seed)
{ // Zilong Tan's fast-hash (2012). Modified for byte lengths and portable, unaligned little-endian reads.
  const uint8_t *bytes = key;
  const uint64_t multiplier = UINT64_C(0x880355f21e6d1965);
  uint64_t hash;
  memcpy (&hash, seed, sizeof (hash));
  hash ^= length * multiplier;
  size_t position = 0;
  for (; length - position >= 8; position += 8)
    hash = (hash ^ crpx_hash_fastmix64 (crpx_hash_read64 (bytes + position))) * multiplier;
  if (position < length) {
    uint64_t tail = 0;
    for (size_t i = 0; position + i < length; i++) tail |= (uint64_t)bytes[position + i] << (8 * i);
    hash = (hash ^ crpx_hash_fastmix64 (tail)) * multiplier;
  }
  return crpx_hash_fastmix64 (hash);
}

uint32_t
crpx_hash_jenkins (const void *key, size_t length)
{ // Bob Jenkins' public-domain one-at-a-time hash.
  return jenkins_bytes (key, length, 0);
}

uint32_t
crpx_hash_mailund_seed32 (const void *key, size_t length, void *seed)
{ // The historical name now calls Jenkins' original one-at-a-time algorithm, with an explicit starting value.
  uint32_t initial;
  memcpy (&initial, seed, sizeof (initial));
  return jenkins_bytes (key, length, initial);
}

uint32_t
crpx_hash_jenkins_mailund_seed32 (const void *key, size_t length, void *seed)
{ // The historical name now calls Bob Jenkins' original lookup2 (1996), replacing the Mailund adaptation.
  uint32_t initial;
  memcpy (&initial, seed, sizeof (initial));
  return jenkins_lookup2 (key, length, initial);
}

uint32_t
crpx_hash_fletcher32 (const void *key, size_t length)
{ // Fletcher-32 sums little-endian 16-bit words modulo 65535; pad an odd final byte with zero.
  const uint8_t *bytes = key;
  uint32_t first = 0, second = 0;
  size_t position = 0;
  while (length - position >= 2) {
    size_t words = (length - position) / 2;
    if (words > 360) words = 360;
    for (size_t i = 0; i < words; i++, position += 2) {
      first += crpx_hash_read16 (bytes + position);
      second += first;
    }
    first %= 65535;
    second %= 65535;
  }
  if (position < length) {
    first = (first + bytes[position]) % 65535;
    second = (second + first) % 65535;
  }
  return second << 16 | first;
}

uint64_t
crpx_hash_pearson_seed2048 (const void *key, size_t length, const void *table)
{ // Eight Pearson recurrences (Pearson, CACM 33:677, 1990); the table contains a permutation of 0..255.
  const uint8_t *bytes = key, *permutation = table;
  if (!length) return 0;
  uint64_t result = 0;
  for (unsigned lane = 0; lane < 8; lane++) {
    uint8_t hash = permutation[(bytes[0] + lane) & 255];
    for (size_t i = 1; i < length; i++) hash = permutation[hash ^ bytes[i]];
    result |= (uint64_t)hash << (8 * lane);
  }
  return result;
}

uint32_t
crpx_hash_pseudocrc32_seed8192 (const void *key, size_t length, const void *table, uint32_t crc)
{ // A table-driven byte recurrence; the caller's table need not describe a CRC polynomial.
  const uint8_t *bytes = key, *entries = table;
  uint32_t hash = ~crc;
  for (size_t i = 0; i < length; i++) {
    uint32_t entry;
    memcpy (&entry, entries + 4 * ((hash ^ bytes[i]) & 255), sizeof (entry));
    hash = entry ^ (hash >> 8);
  }
  return ~hash;
}

/* Copyright (c) 2010 Paul Hsieh. BSD-3-Clause; see LICENSES/SuperFastHash-BSD-3-Clause.txt.
   Modified for Curupixa: size_t byte counts, explicit seed, unsigned tails and portable loads. */
uint32_t
crpx_hsieh_hash32_seed32 (const void *key, size_t len, void *seed)
{
  const uint8_t *data = key;
  uint32_t hash, tmp;
  memcpy (&hash, seed, sizeof (hash));
  size_t rem;

  rem = len & 3;
  len >>= 2;

  /* Main loop */
  for (; len > 0; len--) {
    hash  += crpx_hash_read16 (data);
    tmp    = (uint32_t)(crpx_hash_read16 (data+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    data  += 2*sizeof (uint16_t);
    hash  += hash >> 11;
  }

  /* Handle end cases */
  switch (rem) {
    case 3:
      hash += crpx_hash_read16 (data);
      hash ^= hash << 16;
      hash ^= (uint32_t)data[2] << 18;
      hash += hash >> 11;
      break;
    case 2:
      hash += crpx_hash_read16 (data);
      hash ^= hash << 11;
      hash += hash >> 17;
      break;
    case 1:
      hash += (uint32_t)*data;
      hash ^= hash << 10;
      hash += hash >> 1;
  }

  /* Final avalanche. */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}

uint64_t
crpx_murmurhash3_128bits (const void *key, const size_t len, const uint32_t seed, void *out) // out[] is 128 bits
{ // MurmurHash3 from Austin Appleby's public-domain SMHasher MurmurHash3_x64_128; see docs/random-hashing.md
  const uint8_t *data = key ? key : (const void *)"";
  const size_t nblocks = len / 16;
  size_t i;
  uint64_t h1 = seed, h2 = seed;

  for(i = 0; i < nblocks; i++) {
    uint64_t k1 = crpx_hash_read64 (data + i*16);
    uint64_t k2 = crpx_hash_read64 (data + i*16 + 8);
    k1 *= 0x87c37b91114253d5ULL; k1 = (k1 << 31) | (k1 >> 33); k1 *= 0x4cf5ad432745937fULL; h1 ^= k1;
    h1 = (h1 << 27) | (h1 >> 37); h1 += h2; h1 = h1 * 5 + 0x52dce729ULL;
    k2 *= 0x4cf5ad432745937fULL; k2 = (k2 << 33) | (k2 >> 31); k2 *= 0x87c37b91114253d5ULL; h2 ^= k2;
    h2 = (h2 << 31) | (h2 >> 33); h2 += h1; h2 = h2 * 5 + 0x38495ab5ULL;
  }
  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
  uint64_t k1 = 0, k2 = 0;

  switch(len & 15) {
    case 15: k2 ^= (uint64_t)(tail[14]) << 48; CRPX_attribute_FALLTHROUGH
    case 14: k2 ^= (uint64_t)(tail[13]) << 40; CRPX_attribute_FALLTHROUGH
    case 13: k2 ^= (uint64_t)(tail[12]) << 32; CRPX_attribute_FALLTHROUGH
    case 12: k2 ^= (uint64_t)(tail[11]) << 24; CRPX_attribute_FALLTHROUGH
    case 11: k2 ^= (uint64_t)(tail[10]) << 16; CRPX_attribute_FALLTHROUGH
    case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;  CRPX_attribute_FALLTHROUGH
    case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;
             k2 *= 0x4cf5ad432745937fULL; k2 = (k2 << 33) | (k2 >> 31);
             k2 *= 0x87c37b91114253d5ULL; h2 ^= k2; CRPX_attribute_FALLTHROUGH
    case  8: k1 ^= (uint64_t)(tail[ 7]) << 56; CRPX_attribute_FALLTHROUGH
    case  7: k1 ^= (uint64_t)(tail[ 6]) << 48; CRPX_attribute_FALLTHROUGH
    case  6: k1 ^= (uint64_t)(tail[ 5]) << 40; CRPX_attribute_FALLTHROUGH
    case  5: k1 ^= (uint64_t)(tail[ 4]) << 32; CRPX_attribute_FALLTHROUGH
    case  4: k1 ^= (uint64_t)(tail[ 3]) << 24; CRPX_attribute_FALLTHROUGH
    case  3: k1 ^= (uint64_t)(tail[ 2]) << 16; CRPX_attribute_FALLTHROUGH
    case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;  CRPX_attribute_FALLTHROUGH
    case  1: k1 ^= (uint64_t)(tail[ 0]) << 0; // equiv to "case 1"
             k1 *= 0x87c37b91114253d5ULL; k1 = (k1 << 31) | (k1 >> 33); k1 *= 0x4cf5ad432745937fULL; h1 ^= k1; break;
  };

  h1 ^= len; h2 ^= len;
  h1 += h2; h2 += h1;
  h1 ^= h1 >> 33; h1 *= 0xff51afd7ed558ccdULL; h1 ^= h1 >> 33; h1 *= 0xc4ceb9fe1a85ec53ULL; h1 ^= h1 >> 33;
  h2 ^= h2 >> 33; h2 *= 0xff51afd7ed558ccdULL; h2 ^= h2 >> 33; h2 *= 0xc4ceb9fe1a85ec53ULL; h2 ^= h2 >> 33;
  h1 += h2; h2 += h1;
  if (out) {
    uint64_t result[2] = {h1, h2};
    memcpy (out, result, sizeof (result));
  }
  return crpx_mumhash64_mixer (h1, h2);
}

uint32_t
crpx_murmurhash3_32bits (const void *data, const size_t nbytes, const uint32_t seed)
{ // Austin Appleby's public-domain SMHasher MurmurHash3_x86_32.
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const size_t nblocks = nbytes / 4;
  const uint8_t *bytes = data ? data : (const void *)"";
  const uint8_t *tail = bytes + (nblocks * 4);
  uint32_t k, h = seed;
  size_t i;

  for (i = 0; i < nblocks; i++) {
    k = crpx_hash_read32 (bytes + i*4); k *= c1; k = (k << 15) | (k >> (32 - 15)); k *= c2;
    h ^= k; h = (h << 13) | (h >> (32 - 13)); h = (h * 5) + 0xe6546b64;
  }
  k = 0;
  switch (nbytes & 3) {
    case 3: k ^= tail[2] << 16; CRPX_attribute_FALLTHROUGH
    case 2: k ^= tail[1] << 8;  CRPX_attribute_FALLTHROUGH
    case 1: k ^= tail[0]; k *= c1; k = (k << 15) | (k >> (32 - 15)); k *= c2; h ^= k; break;
  };
  h ^= nbytes; h ^= h >> 16; h *= 0x85ebca6b; h ^= h >> 13; h *= 0xc2b2ae35; h ^= h >> 16;
  return h;
}

/* SipHash reference C implementation:  https://github.com/veorq/SipHash  (Public Domain)
   Copyright (c) 2012-2021 Jean-Philippe Aumasson <jeanphilippe.aumasson@gmail.com>
   Copyright (c) 2012-2014 Daniel J. Bernstein <djb@cr.yp.to> */
#define U8TO64_LE(p)                                                           \
    (((uint64_t)((p)[0])) | ((uint64_t)((p)[1]) << 8) |                        \
     ((uint64_t)((p)[2]) << 16) | ((uint64_t)((p)[3]) << 24) |                 \
     ((uint64_t)((p)[4]) << 32) | ((uint64_t)((p)[5]) << 40) |                 \
     ((uint64_t)((p)[6]) << 48) | ((uint64_t)((p)[7]) << 56))
#define SIPROUND                                                               \
    do {                                                                       \
        v0 += v1; v1 = ROTL64(v1, 13); v1 ^= v0; v0 = ROTL64(v0, 32);          \
        v2 += v3; v3 = ROTL64(v3, 16); v3 ^= v2;                               \
        v0 += v3; v3 = ROTL64(v3, 21); v3 ^= v0;                               \
        v2 += v1; v1 = ROTL64(v1, 17); v1 ^= v2; v2 = ROTL64(v2, 32);          \
    } while (0)

uint64_t
crpx_siphash128_seed128 (const void *in, const size_t inlen, const void *seed, void *out)
{ // k is 16 bytes seed (128bits), out 128bits (use hiphash64 if you're only interested in 64 bits)
  const uint8_t *ni = in ? in : (const void *)"";
  const uint8_t *kk = (const uint8_t*) seed;
  uint64_t v0 = 0x736f6d6570736575ULL, v1 = 0x646f72616e646f6dULL, v2 = 0x6c7967656e657261ULL, v3 = 0x7465646279746573ULL;
  uint64_t k0 = U8TO64_LE(kk);
  uint64_t k1 = U8TO64_LE(kk + 8);
  uint64_t m, result[2];
  uint8_t i, c_rounds = 2, d_rounds = 4; // main constants for algo, this means siphash-2-4
  const uint8_t *end = ni + inlen - (inlen % sizeof(uint64_t));
  const int left = inlen & 7;
  uint64_t b = ((uint64_t)inlen) << 56;
  v3 ^= k1;
  v2 ^= k0;
  v1 ^= k1;
  v0 ^= k0;
  v1 ^= 0xee; // distinct from hiphash64
  for (; ni != end; ni += 8) {
    m = U8TO64_LE(ni);
    v3 ^= m;
    for (i = 0; i < c_rounds; ++i) SIPROUND;
    v0 ^= m;
  }
  switch (left) {
    case 7: b |= ((uint64_t)ni[6]) << 48; CRPX_attribute_FALLTHROUGH
    case 6: b |= ((uint64_t)ni[5]) << 40; CRPX_attribute_FALLTHROUGH
    case 5: b |= ((uint64_t)ni[4]) << 32; CRPX_attribute_FALLTHROUGH
    case 4: b |= ((uint64_t)ni[3]) << 24; CRPX_attribute_FALLTHROUGH
    case 3: b |= ((uint64_t)ni[2]) << 16; CRPX_attribute_FALLTHROUGH
    case 2: b |= ((uint64_t)ni[1]) << 8;  CRPX_attribute_FALLTHROUGH
    case 1: b |= ((uint64_t)ni[0]);  break;
    case 0: break;
  }
  v3 ^= b;
  for (i = 0; i < c_rounds; ++i) SIPROUND;
  v0 ^= b;
  v2 ^= 0xee; // distinct from siphash64
  for (i = 0; i < d_rounds; ++i) SIPROUND;
  result[0] = v0 ^ v1 ^ v2 ^ v3;
  v1 ^= 0xdd;
  for (i = 0; i < d_rounds; ++i) SIPROUND;
  result[1] = v0 ^ v1 ^ v2 ^ v3;
  if (out) memcpy(out, result, 16);
  return crpx_wyhash64_mixer (result[0], result[1]);
}

uint64_t
crpx_siphash64_seed128 (const void *in, const size_t inlen, const void *seed)
{ // k is 16 bytes seed (128bits)
  const uint8_t *ni = in ? in : (const void *)"";
  const uint8_t *kk = (const uint8_t*) seed;
  uint64_t v0 = 0x736f6d6570736575ULL, v1 = 0x646f72616e646f6dULL, v2 = 0x6c7967656e657261ULL, v3 = 0x7465646279746573ULL;
  uint64_t k0 = U8TO64_LE(kk);
  uint64_t k1 = U8TO64_LE(kk + 8);
  uint64_t m;
  uint8_t i, c_rounds = 2, d_rounds = 4; // main constants for algo, this means siphash-2-4
  const uint8_t *end = ni + inlen - (inlen % sizeof(uint64_t));
  const int left = inlen & 7;
  uint64_t b = ((uint64_t)inlen) << 56;
  v3 ^= k1;
  v2 ^= k0;
  v1 ^= k1;
  v0 ^= k0;
  for (; ni != end; ni += 8) {
    m = U8TO64_LE(ni);
    v3 ^= m;
    for (i = 0; i < c_rounds; ++i) SIPROUND;
    v0 ^= m;
  }
  switch (left) {
    case 7: b |= ((uint64_t)ni[6]) << 48; CRPX_attribute_FALLTHROUGH
    case 6: b |= ((uint64_t)ni[5]) << 40; CRPX_attribute_FALLTHROUGH
    case 5: b |= ((uint64_t)ni[4]) << 32; CRPX_attribute_FALLTHROUGH
    case 4: b |= ((uint64_t)ni[3]) << 24; CRPX_attribute_FALLTHROUGH
    case 3: b |= ((uint64_t)ni[2]) << 16; CRPX_attribute_FALLTHROUGH
    case 2: b |= ((uint64_t)ni[1]) << 8;  CRPX_attribute_FALLTHROUGH
    case 1: b |= ((uint64_t)ni[0]);  break;
    case 0: break;
  }
  v3 ^= b;
  for (i = 0; i < c_rounds; ++i) SIPROUND;
  v0 ^= b;
  v2 ^= 0xff;
  for (i = 0; i < d_rounds; ++i) SIPROUND;
  return v0 ^ v1 ^ v2 ^ v3;
}

static uint32_t
jenkins_bytes (const uint8_t *key, size_t length, uint32_t hash)
{
  for (size_t i = 0; i < length; i++) {
    hash += key[i];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  hash += hash << 3;
  hash ^= hash >> 11;
  return hash + (hash << 15);
}

static uint32_t
jenkins_lookup2 (const uint8_t *key, size_t length, uint32_t hash)
{ // Adapted from Bob Jenkins' lookup2.c (December 1996, public domain), including total length before the tail.
  uint32_t a = UINT32_C(0x9e3779b9), b = a, c = hash;
  size_t remaining = length;
#define JENKINS_MIX() do { \
  a -= b; a -= c; a ^= c >> 13; b -= c; b -= a; b ^= a << 8; c -= a; c -= b; c ^= b >> 13; \
  a -= b; a -= c; a ^= c >> 12; b -= c; b -= a; b ^= a << 16; c -= a; c -= b; c ^= b >> 5; \
  a -= b; a -= c; a ^= c >> 3; b -= c; b -= a; b ^= a << 10; c -= a; c -= b; c ^= b >> 15; \
} while (0)
  while (remaining >= 12) {
    a += crpx_hash_read32 (key);
    b += crpx_hash_read32 (key + 4);
    c += crpx_hash_read32 (key + 8);
    JENKINS_MIX ();
    key += 12;
    remaining -= 12;
  }
  c += (uint32_t)length;
  for (size_t i = 0; i < remaining; i++) {
    if (i < 4) a += (uint32_t)key[i] << (8 * i);
    else if (i < 8) b += (uint32_t)key[i] << (8 * (i - 4));
    else c += (uint32_t)key[i] << (8 * (i - 7));
  }
  JENKINS_MIX ();
#undef JENKINS_MIX
  return c;
}
