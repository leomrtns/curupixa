/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html). 
 */


/*! \file hash_functions_generators.c 
 *  \brief simple hash functions.  */

#include "hash_functions_generators.h"
#include "internal_random_constants.h" // not available to the user, only locally

inline uint64_t
crpx_mumhash64_mixer (uint64_t a, uint64_t b)
{ 
  uint64_t ha = a >> 32, hb = b >> 32, la = a & 0xFFFFFFFF00000000ULL, lb = b & 0xFFFFFFFF00000000ULL;
  uint64_t rh = ha * hb, rm_0 = ha * lb, rm_1 = hb * la, rl =  la * lb, t = rl + (rm_0<<32);
  a = t + (rm_1 << 32); b = rh + (rm_0 >> 32) + (rm_1 >> 32); a += b;
  return a;
}

inline uint64_t
crpx_wyhash64_mixer (uint64_t a, uint64_t b)
{  // from  Wang Yi wyhash, very similar to mumhash
  uint64_t ha = a >> 32, hb = b >> 32, la = a & 0xFFFFFFFF00000000ULL, lb = b & 0xFFFFFFFF00000000ULL;
  uint64_t rh = ha * hb, rm_0 = ha * lb, rm_1 = hb * la, rl =  la * lb, t = rl + (rm_0<<32), c = t < rl;
  b ^= 0x60bee2bee120fc15ull;
  a = t + (rm_1 << 32); c += a < t; b = rh + (rm_0 >> 32) + (rm_1 >> 32) + c; a ^= b;
  return a;
}

inline uint32_t 
crpx_hash_64_to_32 (uint64_t key)
{ // https://gist.github.com/badboy/6267743#64-bit-to-32-bit-hash-functions 
  key = (~key) + (key << 18); key = key ^ (key >> 31); key = key * 21;
  key = key ^ (key >> 11);    key = key + (key << 6);  key = key ^ (key >> 22);
  return (uint32_t) key;
}

/* single integer hash functions */

inline uint64_t 
crpx_hashint_staffordmix64 (uint64_t z) 
{ // biomcmc and https://github.com/vigna/MRG32k3a/blob/master/MRG32k3a.c (aka Mix13 by David Stafford see blogpost below)
  z += 0xbd63743fULL; // 32bit primer number added by @leomrtns; o.w. identical to splitmix64
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9; z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

inline uint64_t
crpx_hashint_splitmix64 (uint64_t x) // same as rng_splitmix with state=0; hashint_staffordmix is this with a different state 
{ // https://github.com/skeeto/hash-prospector 
  x ^= x >> 30; x *= 0xbf58476d1ce4e5b9U; x ^= x >> 27; x *= 0x94d049bb133111ebU; 
  return x ^ (x >> 31);
}

inline uint64_t
crpx_hashint_splitmix64_inverse (uint64_t x)
{ // https://github.com/skeeto/hash-prospector 
  x ^= x >> 31 ^ x >> 62; x *= 0x319642b2d24d8ec3ULL;
  x ^= x >> 27 ^ x >> 54; x *= 0x96de1b173f119089ULL; x ^= x >> 30 ^ x >> 60;
  return x;
}

inline uint64_t
crpx_hashint_degski64 (uint64_t x)
{ // https://github.com/skeeto/hash-prospector 
  // x ^= 0x1e018aaf2b12443ULL; // 64 bit prime number added by @leomrtns does not increase dieharder/practrand performance
  x ^= x >> 32; x *= 0xd6e8feb86659fd93ULL;
  x ^= x >> 32; x *= 0xd6e8feb86659fd93ULL; x ^= x >> 32;
  return x;
}

inline uint64_t
crpx_hashint_degski64_inverse (uint64_t x)
{ // https://github.com/skeeto/hash-prospector 
  x ^= x >> 32; x *= 0xcfee444d8b59a89bULL;
  x ^= x >> 32; x *= 0xcfee444d8b59a89bULL; x ^= x >> 32;
  //  return x ^ 0x1e018aaf2b12443ULL;
  return x;
}

inline uint64_t
crpx_hashint_fastmix64 (uint64_t x) // terrible dieharder properties; good for compression _not_ PRNGs
{ // https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c; Compression function for Merkle-Damgard construction
  x ^= x >> 23; x *= 0x2127599bf4325c37ULL; x ^= x >> 47;
  return x;
}

inline uint64_t
crpx_hashint_murmurmix64 (uint64_t k) 
{ // https://github.com/lemire/clhash/blob/master/clhash.c ; mixer for murmurhash finaliser 
  k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33;
  return k;
}

inline uint64_t 
crpx_hashint_rrmixer64 (uint64_t x) 
{ // https://mostlymangling.blogspot.com/2018/07/
  x ^= ROTR64(x, 49) ^ ROTR64(x, 24); x *= 0x9fb21c651e98df25LL;
  x ^= x >> 28; x *= 0x9fb21c651e98df25LL;
  return x ^ x >> 28;
}

inline uint64_t 
crpx_hashint_nasam64 (uint64_t x) 
{ // https://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html
  x ^= 0xb50b2ed9ebf398e9ULL; // 64 bits prime number added by @leomrtns
  x ^= ROTR64(x, 25) ^ ROTR64(x, 47); 
  x *= 0x9E6C63D0676A9A99UL; x ^= x >> 23 ^ x >> 51; 
  x *= 0x9E6D62D06F6A9A9BUL; x ^= x >> 23 ^ x >> 51;
  return x;
}

inline uint64_t
crpx_hashint_pelican64 (uint64_t z)
{ // https://github.com/tommyettinger/sarong/blob/master/src/main/java/sarong/PelicanRNG.java
  z ^= 0x9b25c746f0306ff9ULL; // 64 bits prime number added by @leomrtns
  z = (z ^ (z << 41 | z >> 23) ^ (z << 17 | z >> 47) ^ 0xD1B54A32D192ED03ULL) * 0xAEF17502108EF2D9ULL;
  z = (z ^ z >> 43 ^ z >> 31 ^ z >> 23) * 0xDB4F0B9175AE2165ULL;
  return z ^ z >> 28;
}

inline uint64_t 
crpx_hashint_moremur64 (uint64_t x) 
{ // https://mostlymangling.blogspot.com/2019/12/stronger-better-morer-moremur-better.html
  x ^= x >> 27; x *= 0x3C79AC492BA7B653UL;
  x ^= x >> 33; x *= 0x1C69B3F74AC4AE35UL; x ^= x >> 27;
  return x;
}

inline uint64_t
crpx_hashint_entropy (uint64_t x) 
{ // https://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html (by David Stafford)
  x += 0x9a730fb1ULL; // 32bit primer number added by @leomrtns
  x ^= x >> 31; x *= 0x7fb5d329728ea185;
  x ^= x >> 27; x *= 0x81dadef4bc2dd44d; 
  return x ^ (x >> 33);
}

inline uint32_t
crpx_hashint_jenkins (uint32_t a) // slower than most
{ // full avalanche (https://gist.github.com/badboy/6267743) (http://burtleburtle.net/bob/hash/integer.html) and biomcmc
  a = (a+0x7ed55d16) + (a<<12); a = (a^0xc761c23c) ^ (a>>19); a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);  a = (a+0xfd7046c5) + (a<<3);  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

inline uint32_t
crpx_hashint_jenkins_v2 (uint32_t a) // slower than most
{ // full avalanche (http://burtleburtle.net/bob/hash/integer.html) and biomcmc
  a = (a+0x7fb9b1ee) + (a<<12); a = (a^0xab35dd63) ^ (a>>19); a = (a+0x41ed960d) + (a<<5); 
  a = (a+0xc7d0125e) ^ (a<<9);  a = (a+0x071f9f8f) + (a<<3);  a = (a^0x55ab55b9) ^ (a>>16);
  return a;
}

inline uint32_t
crpx_hashint_avalanche (uint32_t a) // a bit slower than others
{ // https://burtleburtle.net/bob/hash/integer.html and biomcmc
  a ^= 0xb41bf865U; // 32bit primer number added by @leomrtns
  a -= (a<<6); a ^= (a>>17); a -= (a<<9); a ^= (a<<4); a -= (a<<3); a ^= (a<<10); a ^= (a>>15);
  return a;
}

inline uint32_t
crpx_hashint_murmurmix (uint32_t x) 
{ // mixer for murmurhash32 finaliser
  x ^= x >> 16; x *= 0x85ebca6bU;
  x ^= x >> 13; x *= 0xc2b2ae35U;
  x ^= x >> 16;
  return x;
}

inline uint32_t
crpx_hashint_wellons3ple (uint32_t x)
{ // https://github.com/skeeto/hash-prospector  (minimally biased hashed by Christopher Wellons) 
  x++; // avoid hash(0) = 0 
  x ^= x >> 17; x *= 0xed5ad4bbU;
  x ^= x >> 11; x *= 0xac4c1b51U;
  x ^= x >> 15; x *= 0x31848babU;
  x ^= x >> 14;
  return x;
}

inline uint32_t
crpx_hashint_wellons3ple_inverse (uint32_t x) // inverse of crpx_hashint_triple32()
{ // https://github.com/skeeto/hash-prospector 
  x ^= x >> 14 ^ x >> 28; x *= 0x32b21703U;
  x ^= x >> 15 ^ x >> 30; x *= 0x469e0db1U;
  x ^= x >> 11 ^ x >> 22; x *= 0x79a85073U;
  x ^= x >> 17;
  return x-1;
}

inline uint32_t
crpx_hashint_wellons (uint32_t x)
{ // https://github.com/skeeto/hash-prospector 
  x += 0x34f1U; // 16 bit prime number added by @leomrtns
  x ^= x >> 16; x *= 0x7feb352dU;
  x ^= x >> 15; x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

inline uint32_t
crpx_hashint_wellons_inverse (uint32_t x)
{ // https://github.com/skeeto/hash-prospector 
  x ^= x >> 16; x *= 0x43021123;
  x ^= x >> 15 ^ x >> 30; x *= 0x1d69e2a5;
  x ^= x >> 16;
  return x - 0x34f1U;
}

inline uint32_t
crpx_hashint_degski (uint32_t x)
{ // https://github.com/skeeto/hash-prospector  and https://gist.github.com/degski/6e2069d6035ae04d5d6f64981c995ec2
  // x ^= 0xbfcdd6dU; // 32 bit prime number added by @leomrtns to avoid hash(0) = 0 -- does not increase dieharder/pracrand performance
  x ^= x >> 16; x *= 0X45D9F3B;
  x ^= x >> 16; x *= 0X45D9F3B; x ^= x >> 16;
  return x;
}

inline uint32_t
crpx_hashint_degski_inverse (uint32_t x)
{ // https://github.com/skeeto/hash-prospector and https://gist.github.com/degski/6e2069d6035ae04d5d6f64981c995ec2
  x ^= x >> 16; x *= 0X119DE1F3;
  x ^= x >> 16; x *= 0X119DE1F3; x ^= x >> 16;
  //  return x ^ 0xbfcdd6dU; 
  return x;
}


inline uint16_t 
crpx_hashint_2xor_16bits (uint16_t x) // 2-round xorshift-multiply; bias = 0.00859
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 8; x *= 0x88b5U; x ^= x >> 7; x *= 0xdb2dU; x ^= x >> 9;
  return x;
}

inline uint16_t 
crpx_hashint_3xor_16bits (uint16_t x) // 3-round xorshift-multiply; bias = 0.00459
{ // https://github.com/skeeto/hash-prospector
  x ^= x >>  7; x *= 0x2993U; x ^= x >>  5; x *= 0xe877U; 
  x ^= x >>  9; x *= 0x0235U; x ^= x >> 10;
  return x;
}

inline uint16_t 
crpx_hashint_noxor_16bits (uint16_t x) // No multiplication; bias = 0.02384
{ // https://github.com/skeeto/hash-prospector
  x += x << 7; x ^= x >> 8; x += x << 3; x ^= x >> 2; x += x << 4; x ^= x >> 8;
  return x;
}

/* from 8 bit blocks to 64 bit value */

uint64_t
crpx_hash_pearson_seed2048 (const void *vkey, size_t len, const void *vseed) // seed must have >= 256 bytes
{ // https://github.com/maciejczyzewski/retter/blob/master/algorithms/Pearson/pearson.c
  const uint8_t *key = (const uint8_t *) vkey, *seed = (uint8_t *) vseed;
  uint64_t hash = 0;
  unsigned char h = 0;
  size_t i,j;
  for (j = 0; j < 8; j++) {
    h = seed[(key[0] + j) & 0xff ]; // 0xff = 255
    for (i = 1; i < len; i++) h = seed[h ^ key[i]];
    hash ^= h << (j * 8);
  }
  return hash;
}

uint64_t
crpx_fnv_hash64 (const void* vkey, size_t len) 
{ // https://github.com/opencoff/portable-lib/blob/master/src/fnvhash.c, https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
  const uint8_t *p = (const uint8_t*) vkey; // http://www.isthe.com/chongo/tech/comp/fnv/index.html CC0
  uint64_t h = 0xCBF29CE484222325ULL;
  size_t i;
  for (i = 0; i < len; ++i) {
    h += (h << 1) + (h << 4) + (h << 5) + (h << 7) + (h << 8) + (h << 40);
    h ^= p[i];
  }
  return h;
}

/* from 8 bit blocks to 32 bit value */

uint32_t
crpx_fnv_hash32 (const void* vkey, size_t len)
{ // https://github.com/opencoff/portable-lib/blob/master/src/fnvhash.c, https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
  const uint8_t *p = (const uint8_t*) vkey; // http://www.isthe.com/chongo/tech/comp/fnv/index.html CC0
  uint32_t h = 0x811C9DC5;
  size_t i;
  for (i = 0; i < len; ++i) {
    h += (h<<1) + (h<<4) + (h<<7) + (h<<8) + (h<<24);
    h ^= p[i];
  }
  return h;
}

uint32_t
crpx_hash_pseudocrc32_seed8192 (const void *vkey, size_t len, const void *vseed, uint32_t crc) // seed must have >= 1024 bytes (256 x 32bits), crc can be zero, for chaining
{ // modified from https://github.com/maciejczyzewski/retter/blob/master/algorithms/CRC/crc32.c (NOT CRC32)
  const uint8_t *key = (const uint8_t*) vkey;
  uint32_t *seed = (uint32_t *) vseed;
  crc = crc ^ ~0U;
  while (len--) crc = seed[(crc ^ *key++) & 0xFF] ^ (crc >> 8);
  return crc ^ ~0U;
}

uint32_t 
crpx_hash_jenkins_mailund_seed32 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/mailund/hash/blob/master/HashFunctions/source/hash_strings.c
  uint8_t *input = (uint8_t*) vkey;
  uint32_t a, b; a = b = 0x9e3779b9;
  uint32_t c = *(uint32_t*)(vseed);

  while (len >= 12) {
      a += *((uint32_t*)input);
      b += *((uint32_t*)input + 4);
      c += *((uint32_t*)input + 8);
      MIX32(a,b,c);
      input += 12;
      len -= 12;
  }
  c += len;
  switch(len) {
      case 11: c += input[10] << 24; CRPX_attribute_FALLTHROUGH
      case 10: c += input[9] << 16;  CRPX_attribute_FALLTHROUGH
      case 9 : c += input[8] << 8;   CRPX_attribute_FALLTHROUGH
      case 8 : b += input[7] << 24;  CRPX_attribute_FALLTHROUGH
      case 7 : b += input[6] << 16;  CRPX_attribute_FALLTHROUGH
      case 6 : b += input[5] << 8;   CRPX_attribute_FALLTHROUGH
      case 5 : b += input[4];        CRPX_attribute_FALLTHROUGH
      case 4 : a += input[3] << 24;  CRPX_attribute_FALLTHROUGH
      case 3 : a += input[2] << 16;  CRPX_attribute_FALLTHROUGH
      case 2 : a += input[1] << 8;   CRPX_attribute_FALLTHROUGH
      case 1 : a += input[0]; break;
  }
  MIX32(a,b,c);
  return c;
}

uint32_t
crpx_hash_jenkins (const void *vkey, size_t len) // not bad dieharder
{ // https://github.com/maciejczyzewski/retter/tree/master/algorithms/Jenkins
  const uint8_t *key = (const uint8_t *) vkey;
  uint32_t hash, i;
  for(hash = i = 0; i < len; ++i) {
    hash += key[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
}

uint32_t
crpx_hsieh_hash32_seed32 (const void * vkey, size_t len, void *vseed)
{ // https://github.com/opencoff/portable-lib/blob/master/src/hsieh_hash.c
  const uint8_t* data = (const uint8_t*) vkey;
  uint32_t hash = (*(uint32_t*)(vseed)), tmp;
  size_t rem;

  rem = len & 3;
  len >>= 2;
  for (;len > 0; len--) {
    hash  += (*((uint16_t*)(data)) ); // (uint16_t)(*data) is illegal since void doesn't know size of data[0] 
    tmp    = (*((uint16_t*)(data+2)) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    data  += 4; 
    hash  += hash >> 11;
   }
  switch (rem) {  /* last bits */ 
    case 3: hash  += (*((uint16_t*)(data)) ); 
            hash ^= hash << 16;
            hash ^= data[2] << 18;
            hash += hash >> 11;
            break;
    case 2: hash  += (*((uint16_t*)(data)) ); 
            hash ^= hash << 11;
            hash += hash >> 17;
            break;
    case 1: hash += *data;
            hash ^= hash << 10;
            hash += hash >> 1;
   }
  hash ^= hash << 3;  hash += hash >> 5;/* Force "avalanching" of final bits */
  hash ^= hash << 4;  hash += hash >> 17;
  hash ^= hash << 25; hash += hash >> 6;
  return hash;
}

uint32_t
crpx_hash_mailund_seed32 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/mailund/hash/blob/master/HashFunctions/source/hash_strings.c one_at_a_time
  const uint8_t *input = (const uint8_t*) vkey;
  uint32_t hash = *(uint32_t*)(vseed);
  for (size_t i = 0; i < len; i++) {
      hash += input[i];// combine
      hash += (hash << 10); hash ^= (hash >> 6);// mix
  }
  hash += (hash << 3);  // final mix
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
}

uint32_t
crpx_hash_rotating_seed32 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/mailund/hash/blob/master/HashFunctions/source/hash_strings.c
  const uint8_t *input = (const uint8_t*) vkey;
  uint32_t hash = *(uint32_t*)(vseed);
  for (size_t i = 0; i < len; i++) hash += ROTL32(hash, 4) ^ input[i];
  return hash;
}

/* from 16 or 32 bits to 32 bit value */

uint32_t 
crpx_hash_fletcher32 (const void *vkey, size_t len) // assumes vkey has pair number of number of bytes (o.w. last byte is lost) (terrible dieharder)
{ // https://github.com/maciejczyzewski/retter/tree/master/algorithms/Fletcher
/* The Fletcher checksum cannot distinguish between blocks of all 0  bits and blocks of all 1 bits. 
 * For example, if a 16-bit block in the data word changes from 0x0000 to 0xFFFF, the Fletcher-32 checksum remains the same. */
  const uint16_t *key = (const uint16_t *) vkey;
  uint32_t sum1 = 0xffff, sum2 = 0xffff;
  len <<= 1; // len is in bytes, but we assume 16 bits per key element

  while (len) {
    unsigned tlen = len > 359 ? 359 : len; // in https://github.com/opencoff/portable-lib they use 360
    len -= tlen;
    do { sum2 += sum1 += *key++; } while (--tlen);
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
  }
  /* Second reduction step to reduce sums to 16 bits */
  sum1 = (sum1 & 0xffff) + (sum1 >> 16);
  sum2 = (sum2 & 0xffff) + (sum2 >> 16);
  return sum2 << 16 | sum1;
}

uint64_t
crpx_fasthash64_seed64 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c and https://github.com/drobilla/zix/blob/main/src/digest.c
  const uint64_t m = 0x880355f21e6d1965ULL;
  const uint64_t *pos = (const uint64_t *) vkey, *end = pos + (len << 3);
  uint64_t h = (*(uint64_t*)(vseed)) ^ (len * m);
  uint64_t v;

  for (;pos < end; pos++) h = (h ^ crpx_hashint_fastmix64 (*pos)) * m; // main loop

  const uint8_t *pos2 = (const uint8_t *) pos;
  v = 0;
  switch (len & 7) {
    case 7: v ^= (uint64_t) pos2[6] << 48; CRPX_attribute_FALLTHROUGH 
    case 6: v ^= (uint64_t) pos2[5] << 40; CRPX_attribute_FALLTHROUGH 
    case 5: v ^= (uint64_t) pos2[4] << 32; CRPX_attribute_FALLTHROUGH 
    case 4: v ^= (uint64_t) pos2[3] << 24; CRPX_attribute_FALLTHROUGH 
    case 3: v ^= (uint64_t) pos2[2] << 16; CRPX_attribute_FALLTHROUGH 
    case 2: v ^= (uint64_t) pos2[1] << 8;  CRPX_attribute_FALLTHROUGH 
    case 1: v ^= (uint64_t) pos2[0]; h ^= crpx_hashint_fastmix64 (v); h *= m;
  }
  return crpx_hashint_fastmix64 (h);
}

uint64_t
crpx_metrohash64_v1_seed64 (const void *vkey, size_t vlen, const void *seed) // 32bits seed in original, but cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash64.c
  const uint8_t *ptr = (const uint8_t *) vkey;
  static const uint64_t k0 = 0xC83A91E1ULL, k1 = 0x8648DBDBULL, k2 = 0x7BDEC03BULL, k3 = 0x2F5870A5ULL;
  uint64_t len = (uint64_t) vlen, hash = (((*(uint64_t*)(seed)) + k2) * k0) + len; 
  const uint8_t * const end = ptr + len;

  if (len >= 32) {
    uint64_t v[4];
    v[0] = v[1] = v[2] = v[3] = hash;
    do {
      v[0] += *(uint64_t*)(ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += *(uint64_t*)(ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += *(uint64_t*)(ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += *(uint64_t*)(ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= (end - 32));
    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 33) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 33) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 33) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 33) * k0;
    hash += v[0] ^ v[1];
  }
  if ((end - ptr) >= 16) {
    uint64_t v0 = hash + (*(uint64_t*)(ptr) * k0); ptr += 8; v0 = ROTR64(v0,33) * k1;
    uint64_t v1 = hash + (*(uint64_t*)(ptr) * k1); ptr += 8; v1 = ROTR64(v1,33) * k2;
    v0 ^= ROTR64(v0 * k0, 35) + v1;
    v1 ^= ROTR64(v1 * k3, 35) + v0;
    hash += v1;
  }
  if ((end - ptr) >= 8) { hash += *(uint64_t*)(ptr) * k3; ptr += 8; hash ^= ROTR64(hash, 33) * k1; }
  if ((end - ptr) >= 4) { hash += *(uint32_t*)(ptr) * k3; ptr += 4; hash ^= ROTR64(hash, 15) * k1; }
  if ((end - ptr) >= 2) { hash += *(uint16_t*)(ptr) * k3; ptr += 2; hash ^= ROTR64(hash, 13) * k1; }
  if ((end - ptr) >= 1) { hash += *(uint8_t*)(ptr) * k3; hash ^= ROTR64(hash, 25) * k1; }
  hash ^= ROTR64(hash, 33);
  hash *= k0;
  hash ^= ROTR64(hash, 33);
  return hash;
}

uint64_t
crpx_metrohash64_v2_seed64 (const void *vkey, size_t vlen, const void *seed) // 32bits seed in original, but cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash64.c
  const uint8_t *ptr = (const uint8_t *) vkey;
  static const uint64_t k0 = 0xD6D018F5ULL, k1 = 0xA2AA033BULL, k2 = 0x62992FC1ULL, k3 = 0x30BC5B29UL;
  uint64_t len = (uint64_t) vlen, hash = (((*(uint64_t*)(seed)) + k2) * k0) + len; 
  const uint8_t * const end = ptr + len;

  if (len >= 32) {
    uint64_t v[4];
    v[0] = v[1] = v[2] = v[3] = hash;
    do {
      v[0] += *(uint64_t*)(ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += *(uint64_t*)(ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += *(uint64_t*)(ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += *(uint64_t*)(ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= (end - 32));
    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 30) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 30) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 30) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 30) * k0;
    hash += v[0] ^ v[1];
  }
  if ((end - ptr) >= 16) {
    uint64_t v0 = hash + (*(uint64_t*)(ptr) * k2); ptr += 8; v0 = ROTR64(v0,29) * k3;
    uint64_t v1 = hash + (*(uint64_t*)(ptr) * k2); ptr += 8; v1 = ROTR64(v1,29) * k3;
    v0 ^= ROTR64(v0 * k0, 34) + v1;
    v1 ^= ROTR64(v1 * k3, 34) + v0;
    hash += v1;
  }
  if ((end - ptr) >= 8) { hash += *(uint64_t*)(ptr) * k3; ptr += 8; hash ^= ROTR64(hash, 36) * k1; }
  if ((end - ptr) >= 4) { hash += *(uint32_t*)(ptr) * k3; ptr += 4; hash ^= ROTR64(hash, 15) * k1; }
  if ((end - ptr) >= 2) { hash += *(uint16_t*)(ptr) * k3; ptr += 2; hash ^= ROTR64(hash, 15) * k1; }
  if ((end - ptr) >= 1) { hash += *(uint8_t*)(ptr)  * k3;           hash ^= ROTR64(hash, 23) * k1; }
  hash ^= ROTR64(hash, 28);
  hash *= k0;
  hash ^= ROTR64(hash, 29);
 return hash; 
}

uint64_t
crpx_metrohash128_v1_seed64 (const void *vkey, size_t vlen, const void *seed, void *out) // 32bits seed in original, but always cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash128.c
  const uint8_t *ptr = (const uint8_t *) vkey;
  static const uint64_t k0 = 0xC83A91E1ULL, k1 = 0x8648DBDBULL, k2 = 0x7BDEC03BULL, k3 = 0x2F5870A5ULL;
  uint64_t v[4], len = (uint64_t) vlen; // out can be NULL or have at least 16 bytes (2x64 = 128 bits)
  const uint8_t * const end = ptr + len;

  v[0] = (((*(uint64_t*)(seed)) - k0) * k3) + len;
  v[1] = (((*(uint64_t*)(seed)) + k1) * k2) + len;
  v[2] = v[3] = 0ULL;

  if (len >= 32) {        
    v[2] = (((*(uint64_t*)(seed)) + k0) * k2) + len;
    v[3] = (((*(uint64_t*)(seed)) - k1) * k3) + len;
    do {
      v[0] += *(uint64_t*)(ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += *(uint64_t*)(ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += *(uint64_t*)(ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += *(uint64_t*)(ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= (end - 32));
    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 26) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 26) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 26) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 30) * k0;
  }
  if ((end - ptr) >= 16) {
    v[0] += *(uint64_t*)(ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],33) * k3;
    v[1] += *(uint64_t*)(ptr) * k2; ptr += 8; v[1] = ROTR64(v[1],33) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 17) * k1;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 17) * k0;
  }
  if ((end - ptr) >= 8) { v[0] += *(uint64_t*)(ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],33) * k3; v[0] ^= ROTR64((v[0] * k2) + v[1], 20) * k1; }
  if ((end - ptr) >= 4) { v[1] += *(uint32_t*)(ptr) * k2; ptr += 4; v[1] = ROTR64(v[1],33) * k3; v[1] ^= ROTR64((v[1] * k3) + v[0], 18) * k0; }
  if ((end - ptr) >= 2) { v[0] += *(uint16_t*)(ptr) * k2; ptr += 2; v[0] = ROTR64(v[0],33) * k3; v[0] ^= ROTR64((v[0] * k2) + v[1], 24) * k1; }
  if ((end - ptr) >= 1) { v[1] += *(uint8_t*)(ptr)  * k2;           v[1] = ROTR64(v[1],33) * k3; v[1] ^= ROTR64((v[1] * k3) + v[0], 24) * k0; }
  v[0] += ROTR64((v[0] * k0) + v[1], 13);
  v[1] += ROTR64((v[1] * k1) + v[0], 37);
  v[0] += ROTR64((v[0] * k2) + v[1], 13);
  v[1] += ROTR64((v[1] * k3) + v[0], 37);

  if (out) memcpy(out, v, 16);
  return crpx_mumhash64_mixer (v[0], v[1]);
}

uint64_t
crpx_metrohash128_v2_seed64 (const void *vkey, size_t vlen, const void *seed, void *out) // 32bits seed in original, but always cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash128.c
  const uint8_t *ptr = (const uint8_t *) vkey;
  static const uint64_t k0 = 0xD6D018F5ULL, k1 = 0xA2AA033BULL, k2 = 0x62992FC1ULL, k3 = 0x30BC5B29UL;
  uint64_t v[4], len = (uint64_t) vlen; // vout must have at least 16 bytes (128bits)
  const uint8_t * const end = ptr + len;

  v[0] = (((*(uint64_t*)(seed)) - k0) * k3) + len;
  v[1] = (((*(uint64_t*)(seed)) + k1) * k2) + len;
  v[2] = v[3] = 0ULL;

  if (len >= 32) {
    v[2] = (((*(uint64_t*)(seed)) + k0) * k2) + len;
    v[3] = (((*(uint64_t*)(seed)) - k1) * k3) + len;
    do {
      v[0] += *(uint64_t*)(ptr) * k0; ptr += 8; v[0] = ROTR64(v[0],29) + v[2];
      v[1] += *(uint64_t*)(ptr) * k1; ptr += 8; v[1] = ROTR64(v[1],29) + v[3];
      v[2] += *(uint64_t*)(ptr) * k2; ptr += 8; v[2] = ROTR64(v[2],29) + v[0];
      v[3] += *(uint64_t*)(ptr) * k3; ptr += 8; v[3] = ROTR64(v[3],29) + v[1];
    } while (ptr <= (end - 32));
    v[2] ^= ROTR64(((v[0] + v[3]) * k0) + v[1], 33) * k1;
    v[3] ^= ROTR64(((v[1] + v[2]) * k1) + v[0], 33) * k0;
    v[0] ^= ROTR64(((v[0] + v[2]) * k0) + v[3], 33) * k1;
    v[1] ^= ROTR64(((v[1] + v[3]) * k1) + v[2], 33) * k0;
  }
  if ((end - ptr) >= 16) {
    v[0] += *(uint64_t*)(ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],29) * k3;
    v[1] += *(uint64_t*)(ptr) * k2; ptr += 8; v[1] = ROTR64(v[1],29) * k3;
    v[0] ^= ROTR64((v[0] * k2) + v[1], 29) * k1;
    v[1] ^= ROTR64((v[1] * k3) + v[0], 29) * k0;
  }
  if ((end - ptr) >= 8) { v[0] += *(uint64_t*)(ptr) * k2; ptr += 8; v[0] = ROTR64(v[0],29) * k3; v[0] ^= ROTR64((v[0] * k2) + v[1], 29) * k1; }
  if ((end - ptr) >= 4) { v[1] += *(uint32_t*)(ptr) * k2; ptr += 4; v[1] = ROTR64(v[1],29) * k3; v[1] ^= ROTR64((v[1] * k3) + v[0], 25) * k0; }
  if ((end - ptr) >= 2) { v[0] += *(uint16_t*)(ptr) * k2; ptr += 2; v[0] = ROTR64(v[0],29) * k3; v[0] ^= ROTR64((v[0] * k2) + v[1], 30) * k1; }
  if ((end - ptr) >= 1) { v[1] += *(uint8_t*)(ptr)  * k2;           v[1] = ROTR64(v[1],29) * k3; v[1] ^= ROTR64((v[1] * k3) + v[0], 18) * k0; }
  v[0] += ROTR64((v[0] * k0) + v[1], 33);
  v[1] += ROTR64((v[1] * k1) + v[0], 33);
  v[0] += ROTR64((v[0] * k2) + v[1], 33);
  v[1] += ROTR64((v[1] * k3) + v[0], 33);

  if (out) memcpy(out, v, 16);
  return crpx_mumhash64_mixer (v[0], v[1]);
}

uint64_t
crpx_murmurhash3_128bits (const void *key, const size_t len, const uint32_t seed, void *out) // out[] is 128 bits
{ // MurmurHash3 from biomcmc and https://github.com/PeterScott/murmur3/ and written originally by Austin Appleby (public domain) 
  const uint8_t * data = (const uint8_t*) key;
  const int nblocks = len / 16;
  int i;
  uint64_t h1 = seed, h2 = seed;
  const uint64_t * blocks = (const uint64_t *) (data);

  for(i = 0; i < nblocks; i++) {
    uint64_t k1 = blocks[i*2+0];
    uint64_t k2 = blocks[i*2+1];
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
    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
  }
  return crpx_mumhash64_mixer (h1, h2);
}

uint32_t 
crpx_murmurhash3_32bits (const void *data, const size_t nbytes, const uint32_t seed)
{// biomcmc and https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const int nblocks = nbytes / 4;
  const uint32_t *blocks = (const uint32_t *) (data);
  const uint8_t *tail = (const uint8_t *)(data) + (nblocks * 4);
  uint32_t k, h = seed;
  int i;

  for (i = 0; i < nblocks; i++) {
    k = blocks[i]; k *= c1; k = (k << 15) | (k >> (32 - 15)); k *= c2;
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
  const uint8_t *ni = (const uint8_t*) in;
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
  const uint8_t *ni = (const uint8_t*) in;
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
