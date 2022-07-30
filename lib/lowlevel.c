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

const char prt_col_reset[] = "\e[0m";
const char *prt_col[][8]={ // 0-black   1-red   2-grn   3-yel   4-blu   5-mag   6-cyn   7-white
   {"\e[0;30m",  "\e[0;31m",  "\e[0;32m",  "\e[0;33m",  "\e[0;34m",  "\e[0;35m",  "\e[0;36m",  "\e[0;37m"},  // 0 regular text
   {"\e[1;30m",  "\e[1;31m",  "\e[1;32m",  "\e[1;33m",  "\e[1;34m",  "\e[1;35m",  "\e[1;36m",  "\e[1;37m"},  // 1 regular bold text
   {"\e[4;30m",  "\e[4;31m",  "\e[4;32m",  "\e[4;33m",  "\e[4;34m",  "\e[4;35m",  "\e[4;36m",  "\e[4;37m"},  // 2 regular underline text
   {"\e[40m",    "\e[41m",    "\e[42m",    "\e[43m",    "\e[44m",    "\e[45m",    "\e[46m",    "\e[47m"},    // 3 regular background
   {"\e[0;100m", "\e[0;101m", "\e[0;102m", "\e[0;103m", "\e[0;104m", "\e[0;105m", "\e[0;106m", "\e[0;107m"}, // 4 high intensity underground
   {"\e[0;90m",  "\e[0;91m",  "\e[0;92m",  "\e[0;93m",  "\e[0;94m",  "\e[0;95m",  "\e[0;96m",  "\e[0;97m"},  // 5 high intensity text 
   {"\e[1;90m",  "\e[1;91m",  "\e[1;92m",  "\e[1;93m",  "\e[1;94m",  "\e[1;95m",  "\e[1;96m",  "\e[1;97m"}   // 6 bold high intensity text
};


/* error-safe memory allocation functions */
void *
crpx_malloc (size_t size)
{
  void *value = malloc (size);
  if ((value == NULL) && (size > 0)) printf ("biomcmc_malloc error allocating %lu bites", size);
  return value;
}

void *
crpx_realloc (void *ptr, size_t size)
{
  void *value = (void *) realloc ((void *)ptr, size);
  if ((value == NULL) && (size > 0)) printf ( "biomcmc_realloc error on pointer 0x%p of %lu bites\n", ptr, size);
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

const char *msg_level_colours[] = {"\e[0;101m", "\e[1;31m", "\e[1;33m", "\e[1;34m", "\e[1;32m", "\e[1;37m"};
const char *msg_level_names[] = {"  FATAL", "  ERROR", "WARNING", "   INFO", "VERBOSE", "  DEBUG"};

void
crpx_logger_message (uint8_t level, const char *c_file, int c_line, crpx_global_t cglobal, const char *fmt, ...)
{
  if ((level > cglobal->loglevel_stderr) && (level > cglobal->loglevel_file)) return;

  va_list ap;
  time_t t = time (NULL);
  char msg_prefix[32] = {'\0'};
  strftime (msg_prefix, 32, "%T", localtime (&t)); // alternative is "%F %T" where %F is yyyy-mm-dd and %T is hh:mm:ss 
  /* crpx29 refers to the 30th global variable set (in case it's multithreaded for instance); "%-3u" means to left-adjust, but notice that id can be larger than 3 digits */
  if (level <= cglobal->loglevel_stderr) { // colour output to stderr
    fprintf (stderr, "id%-3u %s %s%s%s ", cglobal->id, msg_prefix, msg_level_colours[level], msg_level_names[level], prt_col_reset);
    va_start (ap, fmt); vfprintf (stderr, fmt, ap); va_end (ap); 
    if ((level < CRPX_LOGLEVEL_WARN) || (level > CRPX_LOGLEVEL_VERBOSE)) fprintf (stderr, "  [file %s line %d]\n", c_file, c_line);
    else fprintf (stderr, "\n");
    fflush(stderr);
  }
  if ((level <= cglobal->loglevel_file) && (cglobal->logfile)) { // no colours to log file
    fprintf (cglobal->logfile, "[id%-3u %s %s] ", cglobal->id, msg_prefix, msg_level_names[level]);
    va_start (ap, fmt); vfprintf (cglobal->logfile, fmt, ap); va_end (ap);
    if ((level < CRPX_LOGLEVEL_WARN) || (level > CRPX_LOGLEVEL_VERBOSE)) fprintf (cglobal->logfile, "  [file %s line %d]\n", c_file, c_line);
    else fprintf (cglobal->logfile, "\n");
    fflush(cglobal->logfile);
  }
  return;
}

void
crpx_logger_set_level (crpx_global_t cglobal, uint8_t level)
{
  if (level > CRPX_LOGLEVEL_DEBUG) level = CRPX_LOGLEVEL_DEBUG;
  cglobal->loglevel_stderr = level;
  crpx_logger_info (cglobal, "Screen log level set to %s", msg_level_names[level]);
  return;
}

void
crpx_logger_set_file (crpx_global_t cglobal, const char *filename, uint8_t level)
{
  if (level > CRPX_LOGLEVEL_DEBUG) level = CRPX_LOGLEVEL_DEBUG;
  if (cglobal->logfile) {
    crpx_logger_warning (cglobal, "crpx_logger_set_file: log file already open, closing it and re-opening as %s", filename);
    fclose (cglobal->logfile);
    cglobal->logfile = NULL;
  }
  cglobal->logfile = fopen (filename, "a");
  if (cglobal->logfile == NULL) {
    crpx_logger_error (cglobal, "crpx_logger_set_file: could not open log file %s", filename);
    return;
  }
  if (cglobal->id) crpx_logger_verbose (cglobal, "crpx_logger_set_file: possible multithreaded application, several threads may write to same log file.");
  cglobal->loglevel_file = level;
  crpx_logger_info (cglobal, "crpx_logger_set_file: file %s opened and log will be appended to it at level %s", filename, msg_level_names[level]);
  return;
}
