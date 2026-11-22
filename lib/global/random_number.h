/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: MIT
 */

/*! \file random_number.h
    \brief Draw random values from a context or a separately seeded stream for each concurrent calculation.
 */
#ifndef _global_random_number_h_
#define _global_random_number_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "context_types.h"
#include "random_number_generators.h"
#include "../status.h"

/* A stream belongs to one worker at a time. The same (seed, stream_id) reproduces its sequence, including seed zero.
   Use an explicit rank/worker identifier rather than thread IDs. Distinct identifiers produce distinct starting
   states; streams are not guaranteed disjoint. For disjoint blocks, copy one seeded stream and jump each copy
   successively; each block allows fewer than 2^128 draws. The generator is xoshiro256++ (not cryptographic). */
typedef struct { uint64_t state[4]; } crpx_rng_stream;
void crpx_rng_stream_init (crpx_rng_stream *stream, uint64_t seed, uint64_t stream_id);
void crpx_rng_stream_jump (crpx_rng_stream *stream);
uint64_t crpx_rng_stream_next (crpx_rng_stream *stream);
crpx_status crpx_rng_stream_range (crpx_rng_stream *stream, uint64_t bound, uint64_t *result);
double crpx_rng_stream_double (crpx_rng_stream *stream);

/* Historical generator IDs 0..20 are retained except retired IDs 2 (Moremur), 13 (RRMixer) and 19 (Xorgens). The context default is ID 17,
   xoshiro256++. Context draws are serialised across arbitrary host threads. Their assignment to workers depends on
   scheduling; use worker-owned streams for reproducible parallel work. Seed zero requests entropy. */
crpx_status crpx_set_random_generator_checked (crpx_global_t context, uint8_t generator, uint64_t seed);
crpx_status crpx_random_range_checked (crpx_global_t context, uint64_t bound, uint64_t *result);
void crpx_set_random_generator (crpx_global_t cglob, uint8_t rng_id, uint64_t seed);
extern uint64_t crpx_random_64bits (crpx_global_t cglob);
extern uint32_t crpx_random_32bits (crpx_global_t cglob);
extern uint32_t crpx_random_32bits_extra (crpx_global_t cglob, uint32_t *extra_result);
extern uint64_t crpx_random_range (crpx_global_t cglob, uint64_t n);
extern double crpx_random_double (crpx_global_t cglob); // [0,1)
extern double crpx_random_double_include_one (crpx_global_t cglob); // [0,1]
extern double crpx_random_double_positive (crpx_global_t cglob); // (0,1)
extern double crpx_random_double_positive_include_one (crpx_global_t cglob); // (0,1]
extern double crpx_random_normal (crpx_global_t cglob, double *extra_result); // generates _2_ random numbers, one is stored in extra_result 
extern double crpx_random_normal_fast (crpx_global_t cglob, double *extra_result); // generates _2_ random numbers, one is stored in extra_result 

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */

