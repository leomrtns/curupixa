/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file hash_functions_generators.c
    Integer permutations and multiplication mixers used to construct hashes. */

/* Integer permutations retain their mathematical constants and Curupixa's historical offsets.
   Wellons: Unlicense; Stafford/SplitMix: public domain; Appleby: public domain; Pelican: Apache-2.0.
   See docs/random-hashing.md for the original sources and scope of the provenance review. */
#include "lowlevel.h"
#include "hash_functions_generators.h"
#include "hash_bytes_internal.h"
#include "internal_random_constants.h" // not available to the user, only locally

uint64_t
crpx_mumhash64_mixer (uint64_t a, uint64_t b)
{
  __uint128_t product = (__uint128_t)a * b;
  return (uint64_t)product + (uint64_t)(product >> 64);
}

uint64_t
crpx_wyhash64_mixer (uint64_t a, uint64_t b)
{  // from  Wang Yi wyhash, very similar to mumhash
  __uint128_t product = (__uint128_t)a * b;
  return (uint64_t)product ^ (uint64_t)(product >> 64);
}

uint32_t
crpx_hash_64_to_32 (uint64_t key)
{ // https://gist.github.com/badboy/6267743#64-bit-to-32-bit-hash-functions
  key = (~key) + (key << 18); key = key ^ (key >> 31); key = key * 21;
  key = key ^ (key >> 11);    key = key + (key << 6);  key = key ^ (key >> 22);
  return (uint32_t) key;
}

/* single integer hash functions */

uint64_t
crpx_hashint_staffordmix64 (uint64_t z)
{ // biomcmc and https://github.com/vigna/MRG32k3a/blob/master/MRG32k3a.c (aka Mix13 by David Stafford see blogpost below)
  z += 0xbd63743fULL; // 32bit primer number added by @leomrtns; o.w. identical to splitmix64
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9; z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

uint64_t
crpx_hashint_splitmix64 (uint64_t x) // same as rng_splitmix with state=0; hashint_staffordmix is this with a different state
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 30; x *= 0xbf58476d1ce4e5b9U; x ^= x >> 27; x *= 0x94d049bb133111ebU;
  return x ^ (x >> 31);
}

uint64_t
crpx_hashint_splitmix64_inverse (uint64_t x)
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 31 ^ x >> 62; x *= 0x319642b2d24d8ec3ULL;
  x ^= x >> 27 ^ x >> 54; x *= 0x96de1b173f119089ULL; x ^= x >> 30 ^ x >> 60;
  return x;
}

uint64_t
crpx_hashint_degski64 (uint64_t x)
{ // https://github.com/skeeto/hash-prospector
  // x ^= 0x1e018aaf2b12443ULL; // 64 bit prime number added by @leomrtns does not increase dieharder/practrand performance
  x ^= x >> 32; x *= 0xd6e8feb86659fd93ULL;
  x ^= x >> 32; x *= 0xd6e8feb86659fd93ULL; x ^= x >> 32;
  return x;
}

uint64_t
crpx_hashint_degski64_inverse (uint64_t x)
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 32; x *= 0xcfee444d8b59a89bULL;
  x ^= x >> 32; x *= 0xcfee444d8b59a89bULL; x ^= x >> 32;
  //  return x ^ 0x1e018aaf2b12443ULL;
  return x;
}

uint64_t
crpx_hashint_fastmix64 (uint64_t x) // terrible dieharder properties; good for compression _not_ PRNGs
{ // Zilong Tan, fast-hash (2012), MIT; see LICENSES/FastHash-MIT.txt.
  return crpx_hash_fastmix64 (x);
}

uint64_t
crpx_hashint_murmurmix64 (uint64_t k)
{ // https://github.com/lemire/clhash/blob/master/clhash.c ; mixer for murmurhash finaliser
  k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33;
  return k;
}

uint64_t
crpx_hashint_rrmixer64 (uint64_t x)
{ // https://mostlymangling.blogspot.com/2018/07/
  x ^= ROTR64(x, 49) ^ ROTR64(x, 24); x *= 0x9fb21c651e98df25LL;
  x ^= x >> 28; x *= 0x9fb21c651e98df25LL;
  return x ^ x >> 28;
}

uint64_t
crpx_hashint_nasam64 (uint64_t x)
{ // https://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html
  x ^= 0xb50b2ed9ebf398e9ULL; // 64 bits prime number added by @leomrtns
  x ^= ROTR64(x, 25) ^ ROTR64(x, 47);
  x *= 0x9E6C63D0676A9A99UL; x ^= x >> 23 ^ x >> 51;
  x *= 0x9E6D62D06F6A9A9BUL; x ^= x >> 23 ^ x >> 51;
  return x;
}

uint64_t
crpx_hashint_pelican64 (uint64_t z)
{ // https://github.com/tommyettinger/sarong/blob/master/src/main/java/sarong/PelicanRNG.java
  z ^= 0x9b25c746f0306ff9ULL; // 64 bits prime number added by @leomrtns
  z = (z ^ (z << 41 | z >> 23) ^ (z << 17 | z >> 47) ^ 0xD1B54A32D192ED03ULL) * 0xAEF17502108EF2D9ULL;
  z = (z ^ z >> 43 ^ z >> 31 ^ z >> 23) * 0xDB4F0B9175AE2165ULL;
  return z ^ z >> 28;
}

uint64_t
crpx_hashint_moremur64 (uint64_t x)
{ // https://mostlymangling.blogspot.com/2019/12/stronger-better-morer-moremur-better.html
  x ^= x >> 27; x *= 0x3C79AC492BA7B653UL;
  x ^= x >> 33; x *= 0x1C69B3F74AC4AE35UL; x ^= x >> 27;
  return x;
}

uint64_t
crpx_hashint_entropy (uint64_t x)
{ // https://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html (by David Stafford)
  x += 0x9a730fb1ULL; // 32bit primer number added by @leomrtns
  x ^= x >> 31; x *= 0x7fb5d329728ea185;
  x ^= x >> 27; x *= 0x81dadef4bc2dd44d;
  return x ^ (x >> 33);
}

uint32_t
crpx_hashint_jenkins (uint32_t a) // slower than most
{ // full avalanche (https://gist.github.com/badboy/6267743) (http://burtleburtle.net/bob/hash/integer.html) and biomcmc
  a = (a+0x7ed55d16) + (a<<12); a = (a^0xc761c23c) ^ (a>>19); a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);  a = (a+0xfd7046c5) + (a<<3);  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

uint32_t
crpx_hashint_jenkins_v2 (uint32_t a) // slower than most
{ // full avalanche (http://burtleburtle.net/bob/hash/integer.html) and biomcmc
  a = (a+0x7fb9b1ee) + (a<<12); a = (a^0xab35dd63) ^ (a>>19); a = (a+0x41ed960d) + (a<<5);
  a = (a+0xc7d0125e) ^ (a<<9);  a = (a+0x071f9f8f) + (a<<3);  a = (a^0x55ab55b9) ^ (a>>16);
  return a;
}

uint32_t
crpx_hashint_avalanche (uint32_t a) // a bit slower than others
{ // https://burtleburtle.net/bob/hash/integer.html and biomcmc
  a ^= 0xb41bf865U; // 32bit primer number added by @leomrtns
  a -= (a<<6); a ^= (a>>17); a -= (a<<9); a ^= (a<<4); a -= (a<<3); a ^= (a<<10); a ^= (a>>15);
  return a;
}

uint32_t
crpx_hashint_murmurmix (uint32_t x)
{ // mixer for murmurhash32 finaliser
  x ^= x >> 16; x *= 0x85ebca6bU;
  x ^= x >> 13; x *= 0xc2b2ae35U;
  x ^= x >> 16;
  return x;
}

uint32_t
crpx_hashint_wellons3ple (uint32_t x)
{ // https://github.com/skeeto/hash-prospector  (minimally biased hashed by Christopher Wellons)
  x++; // avoid hash(0) = 0
  x ^= x >> 17; x *= 0xed5ad4bbU;
  x ^= x >> 11; x *= 0xac4c1b51U;
  x ^= x >> 15; x *= 0x31848babU;
  x ^= x >> 14;
  return x;
}

uint32_t
crpx_hashint_wellons3ple_inverse (uint32_t x) // inverse of crpx_hashint_triple32()
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 14 ^ x >> 28; x *= 0x32b21703U;
  x ^= x >> 15 ^ x >> 30; x *= 0x469e0db1U;
  x ^= x >> 11 ^ x >> 22; x *= 0x79a85073U;
  x ^= x >> 17;
  return x-1;
}

uint32_t
crpx_hashint_wellons (uint32_t x)
{ // https://github.com/skeeto/hash-prospector
  x += 0x34f1U; // 16 bit prime number added by @leomrtns
  x ^= x >> 16; x *= 0x7feb352dU;
  x ^= x >> 15; x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

uint32_t
crpx_hashint_wellons_inverse (uint32_t x)
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 16; x *= 0x43021123;
  x ^= x >> 15 ^ x >> 30; x *= 0x1d69e2a5;
  x ^= x >> 16;
  return x - 0x34f1U;
}

uint32_t
crpx_hashint_degski (uint32_t x)
{ // https://github.com/skeeto/hash-prospector  and https://gist.github.com/degski/6e2069d6035ae04d5d6f64981c995ec2
  // x ^= 0xbfcdd6dU; // 32 bit prime number added by @leomrtns to avoid hash(0) = 0 -- does not increase dieharder/pracrand performance
  x ^= x >> 16; x *= 0X45D9F3B;
  x ^= x >> 16; x *= 0X45D9F3B; x ^= x >> 16;
  return x;
}

uint32_t
crpx_hashint_degski_inverse (uint32_t x)
{ // https://github.com/skeeto/hash-prospector and https://gist.github.com/degski/6e2069d6035ae04d5d6f64981c995ec2
  x ^= x >> 16; x *= 0X119DE1F3;
  x ^= x >> 16; x *= 0X119DE1F3; x ^= x >> 16;
  //  return x ^ 0xbfcdd6dU;
  return x;
}


uint16_t
crpx_hashint_2xor_16bits (uint16_t x) // 2-round xorshift-multiply; bias = 0.00859
{ // https://github.com/skeeto/hash-prospector
  x ^= x >> 8; x *= 0x88b5U; x ^= x >> 7; x *= 0xdb2dU; x ^= x >> 9;
  return x;
}

uint16_t
crpx_hashint_3xor_16bits (uint16_t x) // 3-round xorshift-multiply; bias = 0.00459
{ // https://github.com/skeeto/hash-prospector
  x ^= x >>  7; x *= 0x2993U; x ^= x >>  5; x *= 0xe877U;
  x ^= x >>  9; x *= 0x0235U; x ^= x >> 10;
  return x;
}

uint16_t
crpx_hashint_noxor_16bits (uint16_t x) // No multiplication; bias = 0.02384
{ // https://github.com/skeeto/hash-prospector
  x += x << 7; x ^= x >> 8; x += x << 3; x ^= x >> 2; x += x << 4; x ^= x >> 8;
  return x;
}
