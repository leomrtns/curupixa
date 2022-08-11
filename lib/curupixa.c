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

/*! \file curupixa.c 
 *  \brief high-level functions handling the global variables. 
 */

#include "curupixa.h"

void global_init_logger (crpx_global_t cglob, const char *level_string);
void global_init_simd_instructions (crpx_global_t cglob);
void global_init_threads_rng (crpx_global_t cglob, __attribute__((unused)) uint64_t seed);

crpx_global_t
crpx_global_init (__attribute__((unused)) uint64_t seed, uint16_t thread, const char *level_string)
{
  crpx_global_t cglob = (crpx_global_t) malloc (sizeof (crpx_global_struct));
  if (cglob == NULL) {  fprintf (stderr, "toplevel FATAL ERROR: could not allocate memory for crpx_global_t\n");  return NULL; }
  cglob->error = false;
  cglob->id = thread;
  crpx_get_time_128bits (cglob->elapsed_time);

  global_init_logger (cglob, level_string);
  global_init_simd_instructions (cglob);
  global_init_threads_rng (cglob, seed);

  return cglob;
}

void
global_init_logger (crpx_global_t cglob, const char *level_string)
{
  char level_stdout[16] = {'\0'};
  cglob->loglevel_file = CRPX_LOGLEVEL_DEBUG + 1; // debug is highest, so never uses file unless crpx_set_file() is called
  cglob->logfile = NULL;
  switch(level_string[0]) {
    case 'f': case 'F': cglob->loglevel_stderr = CRPX_LOGLEVEL_FATAL; strcpy(level_stdout,"fatal"); break;
    case 'e': case 'E': cglob->loglevel_stderr = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
    case 'w': case 'W': cglob->loglevel_stderr = CRPX_LOGLEVEL_WARN; strcpy(level_stdout,"warning"); break;
    case 'i': case 'I': cglob->loglevel_stderr = CRPX_LOGLEVEL_INFO; strcpy(level_stdout,"info"); break;
    case 'v': case 'V': cglob->loglevel_stderr = CRPX_LOGLEVEL_VERBOSE; strcpy(level_stdout,"verbose"); break;
    case 'd': case 'D': cglob->loglevel_stderr = CRPX_LOGLEVEL_DEBUG; strcpy(level_stdout,"debug"); break;
    default: cglob->loglevel_stderr = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
  }
  crpx_logger_verbose (cglob, "Thread-safe global variable set \"id%u\" initialised with log level = %s", cglob->id, level_stdout);
  return;
}

void
global_init_simd_instructions (crpx_global_t cglob)
{
#ifdef __SSE4_2__
  cglob->sse = (__builtin_cpu_supports ("sse4.2") > 0);
  crpx_logger_verbose (cglob, "Compiled with SSE4.2 instructions, which are %s by host machine", cglob->sse ? "enabled" : "disabled");
#endif
#ifdef __AVX2__
  cglob->avx = (__builtin_cpu_supports ("avx2") > 0);
  crpx_logger_verbose (cglob, "Compiled with AVX2 instructions, which are %s by host machine", cglob->avx ? "enabled" : "disabled");
#endif
#if !defined(__SSE4_2__) && !defined(__AVX2__)
  cglob->sse = cglob->avx = false;
  crpx_logger_verbose (cglob, "Compiled without SSE3 or AVX2 instructions, irrespective of host machine capabilities");
#endif
  return;
}

void
global_init_threads_rng (crpx_global_t cglob, __attribute__((unused)) uint64_t seed)
{
#ifdef _OPENMP
  cglob->n_threads = omp_get_max_threads(); // this is set even if user/program does not use threads
  crpx_logger_verbose (cglob, "Multithreading available, which will lead to %u independent random streams", cglob->n_threads);
#else
  cglob->n_threads = 1;
  crpx_logger_verbose (cglob, "Compiled without multithread support");
#endif
  return;
}

void
crpx_global_finalise (crpx_global_t cglob)
{
  if (cglob) {
    if (cglob->logfile) { fclose (cglob->logfile); cglob->logfile = NULL; cglob->loglevel_file = CRPX_LOGLEVEL_DEBUG + 1; }
    crpx_logger_verbose (cglob, "Finalising global variables, program finished in %lf seconds.", crpx_update_elapsed_time_128bits (cglob->elapsed_time));
    free (cglob);
  }
}

