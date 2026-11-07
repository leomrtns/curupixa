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

#include "lowlevel.h" /* configuration/feature macros must precede system headers */
#include "hash_functions.h"
#include "internal_random_constants.h" // not available to the user, only locally

#ifdef CRPX_OS_WINDOWS
int windows_getentropy (void* buf, size_t n);
#endif

size_t
crpx_generate_bytesized_random_seeds_from_cpu (crpx_global_t cglob, void *seed, size_t seed_size)
{
  (void)cglob;
  if (!seed || !seed_size) return 0;
#ifdef CRPX_OS_WINDOWS
  return windows_getentropy(seed, seed_size) == 0 ? seed_size : 0;
#else
  FILE *source = fopen("/dev/urandom", "rb");
  if (!source) return 0;
  size_t count = fread(seed, 1, seed_size, source);
  fclose(source);
  return count;
#endif
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
		if (success == FALSE) { CryptReleaseContext(ctx, 0); return -1; }
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
  (void)cglob;
  if (!seed) return;
  unsigned char *bytes = seed;
  while (seed_size) {
    initial_seed += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = crpx_hashint_splitmix64(initial_seed);
    size_t count = seed_size < sizeof(value) ? seed_size : sizeof(value);
    /* Explicit byte order makes deterministic seeding independent of host alignment and endianness. */
    for (size_t i = 0; i < count; i++) bytes[i] = (unsigned char)(value >> (8 * i));
    bytes += count;
    seed_size -= count;
  }
}

/* time functions */
void
crpx_get_time_128bits (uint64_t time[2])
{
#if _POSIX_TIMERS
  struct timespec now;
  clock_gettime (CLOCK_MONOTONIC, &now);
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
  seconds = ((double)now[0] - (double)past[0]) + ((double)now[1] - (double)past[1]) / TIMEWARP;
  past[0] = now[0]; past[1] = now[1];
  return seconds;
}
