/* Copyright (C) 2022-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file hash_bytes_internal.h
    Read hash keys as little-endian words, without alignment or aliasing assumptions. */
#ifndef CRPX_HASH_BYTES_INTERNAL_H
#define CRPX_HASH_BYTES_INTERNAL_H
#include <stdint.h>
#include <string.h>

static inline uint16_t
crpx_hash_read16 (const void *memory)
{
  const uint8_t *bytes = memory;
  return (uint16_t)((uint16_t)bytes[0] | (uint16_t)bytes[1] << 8);
}

static inline uint32_t
crpx_hash_read32 (const void *memory)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  uint32_t value;
  memcpy (&value, memory, sizeof (value));
  return value;
#else
  const uint8_t *bytes = memory;
  return (uint32_t)bytes[0] | (uint32_t)bytes[1] << 8 | (uint32_t)bytes[2] << 16 | (uint32_t)bytes[3] << 24;
#endif
}

static inline uint64_t
crpx_hash_read64 (const void *memory)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  uint64_t value;
  memcpy (&value, memory, sizeof (value));
  return value;
#else
  const uint8_t *bytes = memory;
  return (uint64_t)crpx_hash_read32 (bytes) | (uint64_t)crpx_hash_read32 (bytes + 4) << 32;
#endif
}
/* Zilong Tan's MIT-licensed fast-hash mixer; inlined here so byte hashes pay no call per word. */
static inline uint64_t
crpx_hash_fastmix64 (uint64_t value)
{
  value ^= value >> 23;
  value *= UINT64_C(0x2127599bf4325c37);
  return value ^ (value >> 47);
}
#endif
