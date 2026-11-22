/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "global/global_variable.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf (stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

/* Fixtures from the original MetroHash v1/v2 and Tan FastHash sources (seed 42), and SipHash's published key 0..15.
   The input is the byte sequence 0..length-1. Source revisions are recorded in docs/random-hashing.md. */
static const struct {
  size_t length;
  uint64_t metro1, metro2, fast, sip;
} fixtures[] = {
  {0, UINT64_C(0x75ff7bea037dbefc), UINT64_C(0x5a37d33a994a641b), UINT64_C(0x7074b3961043c1ef), UINT64_C(0x726fdb47dd0e0e31)},
  {1, UINT64_C(0xdf3c9630f22fbd80), UINT64_C(0x9d989fae61c34281), UINT64_C(0x97ec1b7e8173e077), UINT64_C(0x74f839c593dc67fd)},
  {7, UINT64_C(0x8d71b9e6296a56cd), UINT64_C(0x6094f3e87cc1db67), UINT64_C(0x325764a86c4d58b1), UINT64_C(0xab0200f58b01d137)},
  {8, UINT64_C(0x739850b7d2a3e6ec), UINT64_C(0x5899b58a1d75ba5e), UINT64_C(0x591b20655aa7b939), UINT64_C(0x93f5f5799a932462)},
  {15, UINT64_C(0x7968def81b724081), UINT64_C(0xaa7dde3b12e0d28b), UINT64_C(0x0adac8c77e042207), UINT64_C(0xa129ca6149be45e5)},
  {16, UINT64_C(0x740bb88f21a02f92), UINT64_C(0x7f44d0a0b95e2ca7), UINT64_C(0x96fdd290adb46a12), UINT64_C(0x3f2acc7f57c29bdb)},
  {31, UINT64_C(0x46de4f5fd26ba14b), UINT64_C(0x517e0d4cedb064e3), UINT64_C(0x0e10c29feb300267), UINT64_C(0x32d892fad841c342)},
  {32, UINT64_C(0x9716b56f245f095b), UINT64_C(0x7976294f2ad3c08b), UINT64_C(0x27dacb3299200ac8), UINT64_C(0x7127512f72f27cce)},
  {63, UINT64_C(0x2dd858aa2bdfe638), UINT64_C(0x9e084a28cebe456d), UINT64_C(0x0d1e3808ecd0861d), UINT64_C(0x958a324ceb064572)},
  {64, UINT64_C(0xb3daf5bb983d8aa8), UINT64_C(0x6dc7030f9c2d5d5e), UINT64_C(0x253758dc5b6f0a6c), UINT64_C(0xacd2c40b8502cad8)},
  {127, UINT64_C(0x59b925405e77edb4), UINT64_C(0xa8bbad66dd408006), UINT64_C(0xbb396cf7ae3d896e), UINT64_C(0x9e7db1ba81dc7e3a)},
  {128, UINT64_C(0xacab35fb6dca4b38), UINT64_C(0xe95902273eb4694f), UINT64_C(0xeecc425b7d0b443c), UINT64_C(0xdeb79e256c8736ae)},
};

static int
test_hash_vectors (void)
{
  unsigned char input[128], seed_bytes[17], *unaligned_seed = seed_bytes + 1;
  unsigned char sip_seed[16];
  uint64_t seed = 42;
  for (size_t i = 0; i < sizeof (input); i++) input[i] = (unsigned char)i;
  for (size_t i = 0; i < sizeof (sip_seed); i++) sip_seed[i] = (unsigned char)i;
  memcpy (unaligned_seed, &seed, sizeof (seed));
  for (size_t i = 0; i < sizeof (fixtures)/sizeof (*fixtures); i++) {
    size_t length = fixtures[i].length;
    CHECK (crpx_metrohash64_v1_seed64 (input, length, unaligned_seed) == fixtures[i].metro1);
    CHECK (crpx_metrohash64_v2_seed64 (input, length, unaligned_seed) == fixtures[i].metro2);
    CHECK (crpx_fasthash64_seed64 (input, length, unaligned_seed) == fixtures[i].fast);
    CHECK (crpx_siphash64_seed128 (input, length, sip_seed) == fixtures[i].sip);
  }
  CHECK (crpx_murmurhash3_32bits ("foo", 3, 0) == UINT32_C(0xf6a5c420));
  unsigned char output[17];
  crpx_murmurhash3_128bits ("foo", 3, 0, output + 1);
  uint64_t words[2];
  memcpy (words, output + 1, sizeof (words));
  CHECK (words[0] == UINT64_C(0xe271865701f54561) && words[1] == UINT64_C(0x7eaf87e42bba7d87));
  CHECK (crpx_fnv_hash32 ("hello", 5) == UINT32_C(0xb6fa7167));
  CHECK (crpx_fnv_hash64 ("hello", 5) == UINT64_C(0x7b495389bdbdd4c7));
  CHECK (crpx_mumhash64_mixer (UINT64_MAX, UINT64_MAX) == UINT64_MAX);
  CHECK (crpx_wyhash64_mixer (UINT64_MAX, UINT64_MAX) == UINT64_MAX);
  CHECK (crpx_mumhash64_mixer (1, 2) == 2 && crpx_wyhash64_mixer (1, 2) == 2);
  CHECK (crpx_fasthash64_seed64 (NULL, 0, &seed) == fixtures[0].fast);
  CHECK (crpx_hash_pearson_seed2048 (NULL, 0, sip_seed) == 0);
  return 0;
}

static int
test_byte_bounds (void)
{ // Exact-size allocations make overreads visible under AddressSanitizer; compare all alignment offsets.
  unsigned char input[260], seed_bytes[17], table[256], crc_table[1025];
  uint64_t seed = UINT64_C(0x123456789abcdef0);
  uint32_t seed32 = UINT32_C(0x12345678);
  memcpy (seed_bytes + 1, &seed, sizeof (seed));
  memcpy (seed_bytes + 9, &seed, sizeof (seed));
  for (size_t i = 0; i < sizeof (input); i++) input[i] = (unsigned char)(i * 137 + 19);
  for (size_t i = 0; i < 256; i++) {
    table[i] = (unsigned char)(i * 137 + 19);
    uint32_t entry = (uint32_t)(i * 1234567);
    memcpy (crc_table + 1 + 4*i, &entry, sizeof (entry));
  }
  for (size_t length = 0; length <= 257; length++) for (size_t offset = 0; offset < 8; offset++) {
    unsigned char *allocation = malloc (length + offset + 1);
    CHECK (allocation);
    unsigned char *key = allocation + offset + 1;
    memcpy (key, input, length);
    CHECK (crpx_fasthash64_seed64 (key, length, seed_bytes + 1) == crpx_fasthash64_seed64 (input, length, &seed));
    CHECK (crpx_hash_jenkins (key, length) == crpx_hash_jenkins (input, length));
    CHECK (crpx_hash_jenkins_mailund_seed32 (key, length, &seed32) ==
           crpx_hash_jenkins_mailund_seed32 (input, length, &seed32));
    CHECK (crpx_hash_mailund_seed32 (key, length, &seed32) == crpx_hash_mailund_seed32 (input, length, &seed32));
    CHECK (crpx_hash_fletcher32 (key, length) == crpx_hash_fletcher32 (input, length));
    CHECK (crpx_hash_pearson_seed2048 (key, length, table) == crpx_hash_pearson_seed2048 (input, length, table));
    CHECK (crpx_hash_pseudocrc32_seed8192 (key, length, crc_table + 1, 123) ==
           crpx_hash_pseudocrc32_seed8192 (input, length, crc_table + 1, 123));
    CHECK (crpx_hsieh_hash32_seed32 (key, length, &seed32) == crpx_hsieh_hash32_seed32 (input, length, &seed32));
    CHECK (crpx_fnv_hash32 (key, length) == crpx_fnv_hash32 (input, length));
    CHECK (crpx_fnv_hash64 (key, length) == crpx_fnv_hash64 (input, length));
    CHECK (crpx_murmurhash3_32bits (key, length, 17) == crpx_murmurhash3_32bits (input, length, 17));
    CHECK (crpx_siphash64_seed128 (key, length, seed_bytes + 1) == crpx_siphash64_seed128 (input, length, seed_bytes + 1));
    uint64_t expected[2], actual[2];
    crpx_murmurhash3_128bits (input, length, 17, expected);
    crpx_murmurhash3_128bits (key, length, 17, actual);
    CHECK (!memcmp (expected, actual, sizeof (actual)));
    crpx_siphash128_seed128 (input, length, seed_bytes + 1, expected);
    crpx_siphash128_seed128 (key, length, seed_bytes + 1, actual);
    CHECK (!memcmp (expected, actual, sizeof (actual)));
    CHECK (crpx_metrohash64_v1_seed64 (key, length, seed_bytes + 1) == crpx_metrohash64_v1_seed64 (input, length, &seed));
    CHECK (crpx_metrohash64_v2_seed64 (key, length, seed_bytes + 1) == crpx_metrohash64_v2_seed64 (input, length, &seed));
    crpx_metrohash128_v1_seed64 (input, length, &seed, expected);
    crpx_metrohash128_v1_seed64 (key, length, seed_bytes + 1, actual);
    CHECK (!memcmp (expected, actual, sizeof (actual)));
    crpx_metrohash128_v2_seed64 (input, length, &seed, expected);
    crpx_metrohash128_v2_seed64 (key, length, seed_bytes + 1, actual);
    CHECK (!memcmp (expected, actual, sizeof (actual)));
    free (allocation);
  }
  /* Fletcher's independently calculated modular sums, including an odd byte and multiple reduction blocks. */
  unsigned char all_ones[2001];
  memset (all_ones, 255, sizeof (all_ones));
  CHECK (crpx_hash_fletcher32 (all_ones, 2000) == 0);
  CHECK (crpx_hash_fletcher32 (all_ones, 2001) == UINT32_C(0x00ff00ff));
  return 0;
}

static int
test_streams (void)
{
  crpx_rng_stream stream = {{1, 2, 3, 4}};
  crpx_rng_stream_jump (&stream);
  const uint64_t jumped[] = {
    UINT64_C(0x8c7a153956b5f3d1), UINT64_C(0x701f1a713401d85e),
    UINT64_C(0x6527f66a65469085), UINT64_C(0x8386b786c4408050)
  };
  CHECK (!memcmp (stream.state, jumped, sizeof (jumped)));
  crpx_global_t context = crpx_global_init (123, NULL), reference = crpx_global_init (456, NULL);
  CHECK (context && reference && context->rng_size == 4);
  CHECK (crpx_set_random_generator_checked (reference, 17, 123) == CRPX_OK);
  for (size_t i = 0; i < 100; i++) CHECK (crpx_random_64bits (context) == crpx_random_64bits (reference));
  crpx_global_finalise (context);
  crpx_global_finalise (reference);
  crpx_rng_stream_init (&stream, 0, 0);
  unsigned counts[7] = {0};
  for (size_t i = 0; i < 70000; i++) {
    uint64_t value;
    CHECK (crpx_rng_stream_range (&stream, 7, &value) == CRPX_OK && value < 7);
    counts[value]++;
    CHECK (crpx_rng_stream_range (&stream, UINT64_MAX, &value) == CRPX_OK && value < UINT64_MAX);
    CHECK (crpx_rng_stream_range (&stream, (UINT64_C(1) << 63) + 1, &value) == CRPX_OK);
    CHECK (value < (UINT64_C(1) << 63) + 1);
    CHECK (crpx_hashint_splitmix64_inverse (crpx_hashint_splitmix64 (i)) == i);
    CHECK (crpx_hashint_wellons_inverse (crpx_hashint_wellons ((uint32_t)i)) == i);
  }
  for (size_t i = 0; i < 7; i++) CHECK (counts[i] > 9500 && counts[i] < 10500);
  return 0;
}

int
main (void)
{
  CHECK (!test_hash_vectors ());
  CHECK (!test_byte_bounds ());
  CHECK (!test_streams ());
  return 0;
}
