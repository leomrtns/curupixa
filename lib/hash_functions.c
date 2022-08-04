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


/*! \file hash_functions.c 
 *  \brief simple hash functions.  */

#include "hash_functions.h"
#include "random_constants.h" // not available to the user, only locally

size_t
crpx_generate_bytesized_random_seeds (crpx_global_t cglob, void *seed, size_t seed_size)
{
  size_t i = 0, j = 0, last = seed_size, first = 0;
  int success = 0; // gententropy() and syscall() return 0 on success and negative in failure (syscall returns errno)
  //uint16_t *seed_16bits = (uint16_t *)seed; // up to __builtin_ia32_rdrand64_step() but it tends to have too many zeroes
  uint32_t *seed_32bits = (uint32_t *)seed;

#ifdef HAVE_RDRND 
  for (i=0, j=0; (i < 2 * seed_size) && (j < (seed_size >> 2)); i++) { // DRNG suggests 10 attempts per rng, we do 2
//    success = __builtin_ia32_rdrand16_step (seed_16bits + j); // 16 bits = 2 bytes thus "j+=2" "j<size-1"
    success = __builtin_ia32_rdrand32_step (seed_32bits + j); // 64 bits = 8 bytes thus "j+=8" "j<size-1"
    j += (success > 0); // avoid mispredicted branches 
  } // it may ends with fewer than seed_size since not always succeed (and seed%2 may be > 0)
  j <<= 2; /* j/=4 */
  crpx_logger_verbose (cglob, "Random seeds produced by CPU crystal entropy (RDRAND): %lu", j);
#endif
  first = j;
  seed_size -= j;
  seed += j;
  for (i=0; (seed_size > 0) && (i < 2); i += (success<0)) { // i increments at every failure (which is a negative success)
    j = seed_size > 256 ? 256 : seed_size; /* maximum buffer size is 256 bytes */
#if (__GLIBC__ > 2 || __GLIBC_MINOR__ > 24)
    success = getentropy (seed, j);
#else
    success = syscall (SYS_getrandom, seed, j, 0);
#endif
    seed += j * (~success & 1);      // avoid mispredicted branches (not performance-critical here though)
    seed_size -= j * (~success & 1); // (success==0) is a success 
  }
  if (last-first) crpx_logger_verbose (cglob, "Random seeds produced by linux random: %lu", last - first);
  return last - seed_size;
}

void 
crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4]) /* UNUSED / pilot  */
{
  FILE *fp;
  uint64_t x1[4]; // unsigned long int
  unsigned long long int x = 0; // courtesy of 32bits machines, they're not the same
  unsigned int i, idx, success;
  fp = fopen ("/dev/urandom", "r");
  if (fp == NULL) crpx_logger_info (cglob, "Could not open /dev/urandom to include it in the random seed\n");
  else {
    if (fread (seed, sizeof (uint64_t), 4, fp) != 1) crpx_logger_info (cglob, "Could not read from /dev/urandom to include it in the random seed\n");
    fclose (fp);
  }
  crpx_get_time_128bits (x1);
  seed[0] ^= x1[0];
  seed[1] ^= x1[1];
  seed[2] ^= x1[0];
  seed[3] ^= x1[1];
#ifdef HAVE_RDRND
  for (i=64, idx=0; idx < 4;i--) { // DRNG suggests 10 cycles 
    success = __builtin_ia32_rdrand64_step (&x);
    x1[idx] = x;
    idx += (success > 0); // avoid mispredicted branches (not performance-critical here though) 
  }
#endif
  x = (uint64_t) getpid() + getppid(); // process ID and parent process ID (can be small numbers)
  for (i=0; i<4;i++) seed[i] ^= (x1[i] + ((uint64_t)(i + 1) * x));

  return;
}

/* time functions */
void
crpx_get_time_128bits (uint64_t time[2])
{
#if _POSIX_TIMERS
  struct timespec now;
  clock_gettime (CLOCK_REALTIME, &now);
  time[1] = now.tv_nsec; // always less than 1billion thus 32bits is enough
#else
  struct timeval now;
  gettimeofday (&now, NULL);
  time[1] = now.tv_usec; // always less than 1million thus 32bits is enough
#endif
  time[0] = now.tv_sec;
  return;
}

#ifdef _POSIX_TIMERS
#define TIMEWARP 1.e9
#else
#define TIMEWARP 1.e6
#endif

double
crpx_update_elapsed_time_128bits (uint64_t past[2])
{
  uint64_t now[2];
  double seconds;
  crpx_get_time_128bits (now);
  if (now[1] < past[1]) seconds = (((double)(past[1] - now[1]) / (double)(TIMEWARP)) - 1. + (double)(now[0] - past[0]));
  else                  seconds = (((double)(now[1] - past[1]) / (double)(TIMEWARP))      + (double)(now[0] - past[0]));
  past[0] = now[0]; past[1] = now[1];
  return seconds;
}

/* hash functions */

inline uint64_t
crpx_fastmix64 (uint64_t x) // terrible dieharder properties
{ // https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c; Compression function for Merkle-Damgard construction
  x ^= x >> 23;
  x *= 0x2127599bf4325c37ULL;
  x ^= x >> 47;
  return x;
}

inline uint64_t
crpx_murmurmix64 (uint64_t k) // not bad dieharer properties
{ // https://github.com/lemire/clhash/blob/master/clhash.c ; mixer for murmurhash
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53ULL;
  k ^= k >> 33;
  return k;
}

inline uint64_t 
crpx_rrmixer64 (uint64_t x) 
{ // https://mostlymangling.blogspot.com/2018/07/
  x ^= ROTR64(x, 49) ^ ROTR64(x, 24);
  x *= 0x9fb21c651e98df25LL;
  x ^= x >> 28;
  x *= 0x9fb21c651e98df25LL;
  return x ^ x >> 28;
}

inline uint64_t 
crpx_moremur64 (uint64_t x) 
{ // https://mostlymangling.blogspot.com/2019/12/stronger-better-morer-moremur-better.html
  x ^= x >> 27;
  x *= 0x3C79AC492BA7B653UL;
  x ^= x >> 33;
  x *= 0x1C69B3F74AC4AE35UL;
  x ^= x >> 27;
  return x;
}

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

uint32_t
crpx_hash_pseudocrc32_seed8192 (const void *vkey, size_t len, const void *vseed, uint32_t crc) // seed must have >= 1024 bytes (256 x 32bits), crc can be zero, for chaining
{ // modified from https://github.com/maciejczyzewski/retter/blob/master/algorithms/CRC/crc32.c (NOT CRC32)
  const uint8_t *key = (const uint8_t *) vkey;
  uint32_t *seed = (uint32_t *) vseed;
  crc = crc ^ ~0U;
  while (len--) crc = seed[(crc ^ *key++) & 0xFF] ^ (crc >> 8);
  return crc ^ ~0U;
}

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
crpx_hash_jenkins_mailund_seed32 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/mailund/hash/blob/master/HashFunctions/source/hash_strings.c
  uint32_t a, b; a = b = 0x9e3779b9;
  uint32_t c = *(uint32_t*)(vseed);
  uint8_t *input = (uint8_t*) vkey;

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
crpx_hash_mailund_seed32 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/mailund/hash/blob/master/HashFunctions/source/hash_strings.c one_at_a_time
  uint32_t hash = *(uint32_t*)(vseed);
  const uint8_t *input = (const uint8_t*) vkey;
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
  uint32_t hash = *(uint32_t*)(vseed);
  const uint8_t *input = (const uint8_t*) vkey;
  for (size_t i = 0; i < len; i++) hash += ROTL32(hash, 4) ^ input[i];
  return hash;
}

uint64_t
crpx_fasthash64_seed64 (const void *vkey, size_t len, void *vseed)
{ // https://github.com/opencoff/portable-lib/blob/master/src/fasthash.c
  const uint64_t    m = 0x880355f21e6d1965ULL;
  const uint64_t *pos = (const uint64_t *) vkey;
  const uint64_t *end = pos + (len << 3);
  const unsigned char *pos2;
  uint64_t h = (*(uint64_t*)(vseed)) ^ (len * m);
  uint64_t v;

  for (;pos < end; pos++) { // 8 bytes at a time
    v  = *pos;
    h ^= crpx_fastmix64 (v);
    h *= m;
  }

  pos2 = (const unsigned char*) pos;
  v = 0;
  switch (len & 7) {
    case 7: v ^= (uint64_t) pos2[6] << 48; CRPX_attribute_FALLTHROUGH 
    case 6: v ^= (uint64_t) pos2[5] << 40; CRPX_attribute_FALLTHROUGH 
    case 5: v ^= (uint64_t) pos2[4] << 32; CRPX_attribute_FALLTHROUGH 
    case 4: v ^= (uint64_t) pos2[3] << 24; CRPX_attribute_FALLTHROUGH 
    case 3: v ^= (uint64_t) pos2[2] << 16; CRPX_attribute_FALLTHROUGH 
    case 2: v ^= (uint64_t) pos2[1] << 8;  CRPX_attribute_FALLTHROUGH 
    case 1: v ^= (uint64_t) pos2[0]; h ^= crpx_fastmix64 (v); h *= m;
  }
  return crpx_fastmix64 (h);
}

uint32_t
crpx_fnv_hash32 (const void* vkey, size_t len)
{ // https://github.com/opencoff/portable-lib/blob/master/src/fnvhash.c
  const unsigned char* p = (const unsigned char*) vkey;
  uint32_t h = 2166136261U;
  size_t i;
  for (i = 0; i < len; ++i) {
    h ^= p[i];
    h *= 16777619UL; /* this is 0x0100_0193 */
  }
  return h;
}

uint64_t
crpx_fnv_hash64 (const void* vkey, size_t len) // terrible dieharder properties
{ // https://github.com/opencoff/portable-lib/blob/master/src/fnvhash.c
  const unsigned char* p = (const unsigned char*) vkey;
  uint64_t h = 14695981039346656037ULL;
  size_t i;
  for (i = 0; i < len; ++i) {
    h ^= p[i];
    h *= 1099511628211ULL;  /* this is 0x1000_0000_1b3 */
  }
  return h;
}


uint32_t
crpx_hsieh_hash32_seed32 (const void * vkey, size_t len, void *vseed)
{ // https://github.com/opencoff/portable-lib/blob/master/src/hsieh_hash.c
  const char* data = (const char*) vkey;
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

uint64_t
crpx_metrohash64_v1_seed64 (const void *vkey, size_t vlen, const void *seed) // 32bits seed in original, but cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash64.c
  static const uint64_t k0 = 0xC83A91E1ULL, k1 = 0x8648DBDBULL, k2 = 0x7BDEC03BULL, k3 = 0x2F5870A5ULL;
  const uint8_t *ptr = (const uint8_t *) vkey;
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
  static const uint64_t k0 = 0xD6D018F5ULL, k1 = 0xA2AA033BULL, k2 = 0x62992FC1ULL, k3 = 0x30BC5B29UL;
  const uint8_t *ptr = (const uint8_t *) vkey;
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

void
crpx_metrohash128_v1_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout) // 32bits seed in original, but always cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash128.c
  static const uint64_t k0 = 0xC83A91E1ULL, k1 = 0x8648DBDBULL, k2 = 0x7BDEC03BULL, k3 = 0x2F5870A5ULL;
  const uint8_t *ptr = (const uint8_t *) vkey;
  uint64_t *v = (uint64_t *)vout, len = (uint64_t) vlen; // vout must have at least 16 bytes
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
}

void
crpx_metrohash128_v2_seed64 (const void *vkey, size_t vlen, const void *seed, void *vout) // 32bits seed in original, but always cast to 64bits
{ // https://github.com/opencoff/portable-lib/blob/master/src/metrohash128.c
  static const uint64_t k0 = 0xD6D018F5ULL, k1 = 0xA2AA033BULL, k2 = 0x62992FC1ULL, k3 = 0x30BC5B29UL;
  const uint8_t *ptr = (const uint8_t *) vkey;
  uint64_t *v = (uint64_t *)vout, len = (uint64_t) vlen; // vout must have at least 16 bytes (128bits)
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
}

/* SipHash reference C implementation:  https://github.com/veorq/SipHash  (Public Domain)
   Copyright (c) 2012-2021 Jean-Philippe Aumasson <jeanphilippe.aumasson@gmail.com>
   Copyright (c) 2012-2014 Daniel J. Bernstein <djb@cr.yp.to> */
#define U32TO8_LE(p, v)                                                        \
    (p)[0] = (uint8_t)((v));                                                   \
    (p)[1] = (uint8_t)((v) >> 8);                                              \
    (p)[2] = (uint8_t)((v) >> 16);                                             \
    (p)[3] = (uint8_t)((v) >> 24);
#define U64TO8_LE(p, v)                                                        \
    U32TO8_LE((p), (uint32_t)((v)));                                           \
    U32TO8_LE((p) + 4, (uint32_t)((v) >> 32));
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

void 
crpx_siphash_seed128 (const void *in, const size_t inlen, const void *seed, uint8_t *out, const size_t outlen)
{ // k is 16 bytes seed (128bits), out is 64bits or 128bits 
  assert((outlen == 8) || (outlen == 16)); // 64 or 128 bits
  const unsigned char *ni = (const unsigned char *)in;
  const unsigned char *kk = (const unsigned char *)seed;
  uint64_t v0 = 0x736f6d6570736575ULL, v1 = 0x646f72616e646f6dULL, v2 = 0x6c7967656e657261ULL, v3 = 0x7465646279746573ULL;
  uint64_t k0 = U8TO64_LE(kk);
  uint64_t k1 = U8TO64_LE(kk + 8);
  uint64_t m;
  uint8_t i, c_rounds = 2, d_rounds = 4; // main constants for algo, this means siphash-2-4
  const unsigned char *end = ni + inlen - (inlen % sizeof(uint64_t));
  const int left = inlen & 7;
  uint64_t b = ((uint64_t)inlen) << 56;
  v3 ^= k1;
  v2 ^= k0;
  v1 ^= k1;
  v0 ^= k0;
  if (outlen == 16) v1 ^= 0xee;

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
  if (outlen == 16) v2 ^= 0xee;
  else              v2 ^= 0xff;
  for (i = 0; i < d_rounds; ++i) SIPROUND;
  b = v0 ^ v1 ^ v2 ^ v3;
  U64TO8_LE(out, b);
  if (outlen == 8) return;
  v1 ^= 0xdd;
  for (i = 0; i < d_rounds; ++i) SIPROUND;
  b = v0 ^ v1 ^ v2 ^ v3;
  U64TO8_LE(out + 8, b);
  return;
}
