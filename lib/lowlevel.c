/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.c 
 *  \brief Lowest level basic functions, that should be available to all other modules. 
 */

#include "lowlevel.h"

const char prt_col_reset[] = "\033[0m";
const char *prt_col[][8]={ // 0-black   1-red   2-grn   3-yel   4-blu   5-mag   6-cyn   7-white 
   {"\033[0;30m",  "\033[0;31m",  "\033[0;32m",  "\033[0;33m",  "\033[0;34m",  "\033[0;35m",  "\033[0;36m",  "\033[0;37m"},  // 0 regular text
   {"\033[1;30m",  "\033[1;31m",  "\033[1;32m",  "\033[1;33m",  "\033[1;34m",  "\033[1;35m",  "\033[1;36m",  "\033[1;37m"},  // 1 regular bold text
   {"\033[4;30m",  "\033[4;31m",  "\033[4;32m",  "\033[4;33m",  "\033[4;34m",  "\033[4;35m",  "\033[4;36m",  "\033[4;37m"},  // 2 regular underline text
   {"\033[40m",    "\033[41m",    "\033[42m",    "\033[43m",    "\033[44m",    "\033[45m",    "\033[46m",    "\033[47m"},    // 3 regular background
   {"\033[0;100m", "\033[0;101m", "\033[0;102m", "\033[0;103m", "\033[0;104m", "\033[0;105m", "\033[0;106m", "\033[0;107m"}, // 4 high intensity underground
   {"\033[0;90m",  "\033[0;91m",  "\033[0;92m",  "\033[0;93m",  "\033[0;94m",  "\033[0;95m",  "\033[0;96m",  "\033[0;97m"},  // 5 high intensity text 
   {"\033[1;90m",  "\033[1;91m",  "\033[1;92m",  "\033[1;93m",  "\033[1;94m",  "\033[1;95m",  "\033[1;96m",  "\033[1;97m"}   // 6 bold high intensity text
}; // OBS: '\e' is GNU only; '\033' is pedantic


void // used through a macro (without error message)
crpx_free_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, void *ptr)
{
  if (ptr) { 
    free (ptr); ptr = NULL; 
  } else {
    crpx_logger_message (CRPX_LOGLEVEL_DEBUG, c_file, c_line, cglobal, "Trying to free a NULL pointer.");
  }
  return;
}

void *
crpx_malloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, size_t size)
{
  if (!size) return NULL;
  void *value = malloc (size);
  if (value == NULL) crpx_logger_message (CRPX_LOGLEVEL_ERROR, c_file, c_line, cglobal, 
                                          "CRPX failed memory allocation of %lu bytes: %s. ", size, strerror (errno));
  return value;
}

void *
crpx_calloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, size_t nmemb, size_t size)
{
  if (!nmemb || !size) return NULL;
  void *value = calloc (nmemb, size);
  if (value == NULL) crpx_logger_message (CRPX_LOGLEVEL_ERROR, c_file, c_line, cglobal, 
                                          "CRPX failed memory allocation of %lu members of %lu bytes: %s. ", nmemb, size, strerror (errno));
  return value;
}
// OBS: realloc() and reallocarray() do not touch ptr in case of error. 

void *
crpx_realloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, void *ptr, size_t size)
{
  if (!size) { 
    crpx_logger_message (CRPX_LOGLEVEL_WARN, c_file, c_line, cglobal, "Memory realloc of zero bytes (equivalent to free()) requested at line %d file %s.", c_line, c_file);
    if (ptr) free (ptr); // default behaviour of calloc(ptr, 0) is to free the pointer
    return NULL;
  }
  void *value = (void *) realloc ((void *)ptr, size); // will free ptr if pointer area has moved and return pointer to new area
  if (value == NULL) { // realloc() behaviour is to return NULL but keep ptr intact (in case of overflow user may solve by creating several vectors for instance)
    crpx_logger_message (CRPX_LOGLEVEL_ERROR, c_file, c_line, cglobal, 
                         "CRPX failed memory reallocation of %lu bytes: %s. ", size, strerror (errno));
    // if (ptr) free (ptr); // alternative solution for realloc() and reallocarray(), but prevents user from realising that it failed
  }
  return value;
}

void *
crpx_reallocarray_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, void *ptr, size_t nmemb, size_t size)
{
  if (!size || !nmemb) {
      crpx_logger_message (CRPX_LOGLEVEL_WARN, c_file, c_line, cglobal, // c_file and c_line are _not_ printed by logger in WARNING level
                           "Memory reallocarray of %lu members of %lu bytes (equivalent to free()) requested at line %d file %s.", nmemb, size, c_line, c_file);
    if (ptr) free (ptr); // default behaviour of calloc(ptr, 0) is to free the pointer
    return NULL;
  }
#ifndef CRPX_OS_WINDOWS
  void *value = (void *) reallocarray ((void *)ptr, nmemb, size); // will free ptr if pointer area has moved and return pointer to new area
#else
  #define MUL_NO_OVERFLOW ((size_t)1 << (sizeof(size_t) * 4))
  void *value = NULL;
  if ((nmemb >= MUL_NO_OVERFLOW || size >= MUL_NO_OVERFLOW) && (UINTPTR_MAX / nmemb) < size) { errno = ENOMEM; } 
  else {  value = (void *) realloc ((void *)ptr, nmemb * size); }
#endif  // windows
  if (value == NULL) {
    crpx_logger_message (CRPX_LOGLEVEL_ERROR, c_file, c_line, cglobal, 
                         "CRPX failed memory reallocation_array of %lu members of %lu bytes: %s. ", nmemb, size, strerror (errno));
  }
  return value;
}


void
crpx_fprintf_colour (FILE *stream, int regular, int colour, const char *message, const char *normaltext, ...)
{
  va_list ap;
  if ((regular < 0) || (regular > 6)) regular = 0;
  if ((colour < 0) || (colour > 7)) colour = 1;
  fprintf (stream, "%s%s%s", prt_col[regular][colour], message, prt_col_reset);
  va_start (ap, normaltext);
  vfprintf (stream, normaltext, ap);
  va_end (ap);
}

const char *msg_level_colours[] = {"\033[0;101m", "\033[1;31m", "\033[1;33m", "\033[1;34m", "\033[1;32m", "\033[1;37m"};
const char *msg_level_names[] = {"  FATAL", "  ERROR", "WARNING", "   INFO", "VERBOSE", "  DEBUG"};

void
crpx_logger_message (uint8_t level, const char *c_file, const int c_line, crpx_global_t cglobal, const char *fmt, ...)
{
  if ((level > cglobal->loglevel_stderr) && (level > cglobal->loglevel_file)) return;

  va_list ap;
  uint16_t tid = CRPX_THREAD_NUM; 
  time_t t = time (NULL);
  char msg_prefix[32] = {'\0'};
  #pragma omp critical (crpx_logger_message)
   {
#ifdef CRPX_OS_WINDOWS
    strftime (msg_prefix, 32, "%H:%M:%S", localtime (&t)); // mingw doesn't accept %T but it's identical  
#else
    strftime (msg_prefix, 32, "%T", localtime (&t)); // alternative is "%F %T" where %F is yyyy-mm-dd and %T is hh:mm:ss 
#endif
    /*  "%-3u" means to left-adjust, but notice that thread id can be larger than 3 digits */
    if (level <= cglobal->loglevel_stderr) { // colour output to stderr
      fprintf (stderr, "tid%-3u %s %s%s%s ", tid, msg_prefix, msg_level_colours[level], msg_level_names[level], prt_col_reset);
      va_start (ap, fmt); vfprintf (stderr, fmt, ap); va_end (ap); 
      if ((level < CRPX_LOGLEVEL_WARN) || (level > CRPX_LOGLEVEL_VERBOSE)) fprintf (stderr, "  [file %s line %d]\n", c_file, c_line);
      else fprintf (stderr, "\n");
      fflush(stderr);
    }
    if ((level <= cglobal->loglevel_file) && (cglobal->logfile)) { // no colours to log file
      fprintf (cglobal->logfile, "[tid%-3u %s %s] ", tid, msg_prefix, msg_level_names[level]);
      va_start (ap, fmt); vfprintf (cglobal->logfile, fmt, ap); va_end (ap);
      if ((level < CRPX_LOGLEVEL_WARN) || (level > CRPX_LOGLEVEL_VERBOSE)) fprintf (cglobal->logfile, "  [file %s line %d]\n", c_file, c_line);
      else fprintf (cglobal->logfile, "\n");
      fflush(cglobal->logfile);
    }

    if (level == CRPX_LOGLEVEL_FATAL) cglobal->error = 2; // normal = 0; error = 1; fatal = 2
    else if (level == CRPX_LOGLEVEL_ERROR) cglobal->error = (cglobal->error) ? cglobal->error: 1;
   } // pragma omp critical 
  return;
}

void
crpx_logger_set_level (crpx_global_t cglobal, uint8_t level)
{
  if (level > CRPX_LOGLEVEL_DEBUG) level = CRPX_LOGLEVEL_DEBUG;
  #pragma omp critical (logger_set_level) 
  {
    cglobal->loglevel_stderr = level; 
  }
  crpx_logger_info (cglobal, "Screen log level set to %s", msg_level_names[level]);
  return;
}

void
crpx_logger_set_file (crpx_global_t cglobal, const char *filename, const char *level_string)
{
  char level_stdout[16] = {'\0'};
  #pragma omp single 
   {
    if (cglobal->logfile) {
      crpx_logger_warning (cglobal, "crpx_logger_set_file: log file already open, closing it and re-opening as %s", filename);
      fclose (cglobal->logfile);
      cglobal->logfile = NULL;
    }
    cglobal->logfile = fopen (filename, "a");
    if (cglobal->logfile) {
      cglobal->loglevel_file = crpx_get_logger_level_number (level_string, level_stdout);
      crpx_logger_info (cglobal, "crpx_logger_set_file: file %s opened and log will be appended to it at level %s", filename, level_stdout);
    } else {
      crpx_logger_error (cglobal, "crpx_logger_set_file: could not open log file %s", filename);
    }
   } // pragma omp single 
  return;
}

uint8_t
crpx_get_logger_level_number (const char *level_string, char *level_stdout)
{
  uint8_t loglevel = CRPX_LOGLEVEL_ERROR;
  switch(level_string[0]) {
    case 'f': case 'F': loglevel = CRPX_LOGLEVEL_FATAL; strcpy(level_stdout,"fatal"); break;
    case 'e': case 'E': loglevel = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
    case 'w': case 'W': loglevel = CRPX_LOGLEVEL_WARN; strcpy(level_stdout,"warning"); break;
    case 'i': case 'I': loglevel = CRPX_LOGLEVEL_INFO; strcpy(level_stdout,"info"); break;
    case 'v': case 'V': loglevel = CRPX_LOGLEVEL_VERBOSE; strcpy(level_stdout,"verbose"); break;
    case 'd': case 'D': loglevel = CRPX_LOGLEVEL_DEBUG; strcpy(level_stdout,"debug"); break;
    default: loglevel = CRPX_LOGLEVEL_ERROR; strcpy(level_stdout,"error"); break;
  }
  return loglevel;
}
