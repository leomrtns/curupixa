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
 *  \brief hash functions, higher level in the sense that use curupixa_global_t */

#include "hash_functions.h"
#include "internal_random_constants.h" // not available to the user, only locally

#ifdef CRPX_OS_WINDOWS
int windows_getentropy (void* buf, size_t n);
#endif

size_t
crpx_generate_bytesized_random_seeds_from_cpu (crpx_global_t cglob, void *seed, size_t seed_size)
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
  j <<= 2; /* j*=4 */
  crpx_logger_verbose (cglob, "Number of random bytes produced by CPU crystal entropy (RDRAND): %lu", j);
#endif
  first = j;
  seed_size -= j;
  seed = ((size_t*)seed) + j;
  for (i=0; (seed_size > 0) && (i < 2); i += (success<0)) { // i increments at every failure (which is a negative success)
    j = seed_size > 256 ? 256 : seed_size; /* maximum buffer size is 256 bytes */
#if (__GLIBC__ > 2 || __GLIBC_MINOR__ > 24)
    success = getentropy (seed, j);
#elif defined (CRPX_OS_LINUX) || defined (CRPX_OS_MACOS)
    success = syscall (SYS_getrandom, seed, j, 0);
#else // windows
    success = windows_getentropy (seed, j);
#endif 
    seed = (size_t*)(seed) + (j * (~success & 1));      // avoid mispredicted branches (not performance-critical here though)
    seed_size -= j * (~success & 1); // (success==0) is a success 
  }
  if (last-first) crpx_logger_verbose (cglob, "Number of random bytes produced by OS entropy: %lu", last - first);
  return last - seed_size;
}

#ifdef CRPX_OS_WINDOWS
int 
windows_getentropy (void* buf, size_t n)
{ // https://github.com/dsprenkels/randombytes
	HCRYPTPROV ctx;
	BOOL success;
	DWORD to_read = 0;
	const size_t MAX_DWORD = 0xFFFFFFFF;
	success = CryptAcquireContext(&ctx, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT);
	if (success == FALSE) return -1;

	while (n > 0) {
		to_read = (DWORD)(n < MAX_DWORD ? n : MAX_DWORD);
		success = CryptGenRandom(ctx, to_read, (BYTE*) buf);
		if (success == FALSE) return -1;
		buf = ((char*)buf) + to_read;
		n -= to_read;
	}
	success = CryptReleaseContext(ctx, 0);
	if (success == FALSE) return -1;
	return 0;
}
#endif

void
crpx_generate_bytesized_random_seeds_from_seed (crpx_global_t cglob, void *seed, size_t seed_size, uint64_t initial_seed)
{
  size_t len = seed_size >> 3, rem = seed_size & 7; // len = seed_size%8, rem is remainder
  uint64_t *seed64 = (uint64_t *)seed;
  if (!initial_seed) initial_seed = crpx_list_of_128_random64[ CRPX_THREAD_NUM & 127 ]; // each thread will have its seed
  crpx_logger_verbose (cglob, "Number of random bytes to be produced by seed %lu: %lu", initial_seed, seed_size);
  for (size_t i = 0; i < len; i++) {
    initial_seed += UINT64_C(0x171924dc8e5); // 48bit prime number added by @leomrtns
    seed64[i] = initial_seed = crpx_hashint_moremur64 (initial_seed);
  }
  if (!rem) return;

  initial_seed  = crpx_hashint_moremur64 (initial_seed);
  uint8_t *seed8 = (uint8_t *)(seed);
  seed8 += len << 3;
  switch (rem & 7) {
    case 7: seed8[6] = (uint8_t) (initial_seed >> 48); CRPX_attribute_FALLTHROUGH 
    case 6: seed8[5] = (uint8_t) (initial_seed >> 40); CRPX_attribute_FALLTHROUGH
    case 5: seed8[4] = (uint8_t) (initial_seed >> 32); CRPX_attribute_FALLTHROUGH
    case 4: seed8[3] = (uint8_t) (initial_seed >> 24); CRPX_attribute_FALLTHROUGH
    case 3: seed8[2] = (uint8_t) (initial_seed >> 16); CRPX_attribute_FALLTHROUGH
    case 2: seed8[1] = (uint8_t) (initial_seed >> 8);  CRPX_attribute_FALLTHROUGH
    case 1: seed8[0] = (uint8_t) initial_seed;
  }
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

// for Windows I would need GetSystemTimeAsFileTime() (https://github.com/tezc/sc/blob/master/time/sc_time.c) but I'm assuming mingw64
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
