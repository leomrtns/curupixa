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

crpx_global_t
crpx_global_init (__attribute__((unused)) uint64_t seed, uint16_t thread, const char *level_string)
{
  char level_stdout[16] = {'\0'};
  crpx_global_t cglob = (crpx_global_t) malloc (sizeof (crpx_global_struct));
  if (cglob == NULL) {  fprintf (stderr, "toplevel FATAL ERROR: could not allocate memory for crpx_global_t\n");  return NULL; }
  cglob->loglevel_file = CRPX_LOGLEVEL_DEBUG + 1; // never
  cglob->logfile = NULL;
  cglob->error = false;
  cglob->id = thread;
  switch(level_string[0]) {
    case 'f': case 'F': cglob->loglevel_stderr = CRPX_LOGLEVEL_FATAL; strcpy(level_stdout,"fatal"); break;
    case 'e': case 'E': cglob->loglevel_stderr = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
    case 'w': case 'W': cglob->loglevel_stderr = CRPX_LOGLEVEL_WARN; strcpy(level_stdout,"warning"); break;
    case 'i': case 'I': cglob->loglevel_stderr = CRPX_LOGLEVEL_INFO; strcpy(level_stdout,"info"); break;
    case 'v': case 'V': cglob->loglevel_stderr = CRPX_LOGLEVEL_VERBOSE; strcpy(level_stdout,"verbose"); break;
    case 'd': case 'D': cglob->loglevel_stderr = CRPX_LOGLEVEL_DEBUG; strcpy(level_stdout,"debug"); break;
    default: cglob->loglevel_stderr = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
  }
  crpx_logger_verbose (cglob, "Set id%u of thread-safe global variables initialised with seed=%lu and loglevel=%s", cglob->id, seed, level_stdout);
  return cglob;
}

void
crpx_global_finalise (crpx_global_t cglob)
{
  if (cglob) {
    if (cglob->logfile) { fclose (cglob->logfile); cglob->logfile = NULL; cglob->loglevel_file = CRPX_LOGLEVEL_DEBUG + 1; }
    crpx_logger_verbose (cglob, "Finalising global variables\n");
    free (cglob);
  }
}

