/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
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

/*! \file global_variable.c 
 *  \brief high-level functions handling the global variables. 
 */

#include "global_variable.h"
#include "runtime_internal.h"

crpx_global_t
crpx_global_init (uint64_t seed, const char *level_string)
{
  crpx_global_t cglob = calloc(1, sizeof(*cglob));
  if (!cglob) return NULL;
  cglob->runtime_state = malloc(sizeof(*cglob->runtime_state));
  if (!cglob->runtime_state) { free(cglob); return NULL; }
  atomic_flag_clear(&cglob->runtime_state->rng_lock);
  atomic_flag_clear(&cglob->runtime_state->logger_lock);
  atomic_init(&cglob->runtime_state->references, 1);
  cglob->nthreads = 1;
  cglob->loglevel_file = CRPX_LOGLEVEL_ERROR;
  cglob->loglevel_stderr = crpx_get_logger_level_number(level_string, NULL);
#if defined(__SSE4_2__) && (defined(__GNUC__) || defined(__clang__))
  cglob->sse = (__builtin_cpu_supports("sse4.2") > 0);
#endif
#if defined(__AVX2__) && (defined(__GNUC__) || defined(__clang__))
  cglob->avx = (__builtin_cpu_supports("avx2") > 0);
#endif
  crpx_get_time_128bits(cglob->elapsed_time);
  crpx_set_random_generator(cglob, 0, seed);
  if (!cglob->rng_seed_vector || !cglob->rng_get) {
    crpx_global_finalise(cglob);
    return NULL;
  }
  return cglob;
}

crpx_global_t
crpx_global_retain (crpx_global_t cglob)
{
  if (cglob) atomic_fetch_add_explicit(&cglob->runtime_state->references, 1, memory_order_relaxed);
  return cglob;
}

void
crpx_global_finalise (crpx_global_t cglob)
{
  if (!cglob) return;
  if (atomic_fetch_sub_explicit(&cglob->runtime_state->references, 1, memory_order_acq_rel) != 1) return;
  if (cglob->logfile) fclose(cglob->logfile);
  free(cglob->rng_seed_vector);
  free(cglob->runtime_state);
  free(cglob);
}
