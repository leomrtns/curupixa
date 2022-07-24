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

#include "global_logger.h"

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
biomcmc_malloc (size_t size)
{
  void *value = malloc (size);
  if ((value == NULL) && (size > 0)) biomcmc_error ( "biomcmc_malloc error allocating %d bites", size);
  return value;
}

void *
biomcmc_realloc (void *ptr, size_t size)
{
  void *value = (void *) realloc ((void *)ptr, size);
  if ((value == NULL) && (size > 0)) biomcmc_error ( "biomcmc_realloc error on pointer 0x%08X of %d bites\n", ptr, size);
  return value;
}

void
biomcmc_error (const char *template, ...)
{
  va_list ap;

  //fprintf (stderr, "[%s error] ", PACKAGE_STRING);
  fprintf (stderr, "%s[ error ]%s ", "\e[1;31m", "\e[0m");
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
  fprintf (stderr, "\n");
  fprintf (stderr, "[note to developers] If you want to debug me, set a breakpoint on function biomcmc_error()\n");
  fflush (stderr);
  exit (EXIT_FAILURE);
}

void
biomcmc_warning (const char *template, ...)
{
  va_list ap;

  fprintf (stderr, "%s[warning]%s ", "\e[0;31m", "\e[0m");
  va_start (ap, template);
  vfprintf (stderr, template, ap);
  va_end (ap);
  fprintf (stderr, "\n");
  return;
}

void
biomcmc_fprintf_colour (FILE *stream, int regular, int colour, const char *message, const char *normaltext, ...)
{
  va_list ap;
  if ((regular < 0) || (regular > 6)) regular = 0;
  if ((colour < 0) || (colour > 7)) colour = 1;
  fprintf (stream, "%s%s%s", prt_col[regular][colour], message, prt_col_reset);
  va_start (ap, normaltext);
  vfprintf (stream, normaltext, ap);
  va_end (ap);
}

// DUMP 

static char* getDateString() { // https://github.com/ntpeters/SimpleLogger/blob/master/simplog.c
    // Initialize and get current time
    time_t t = time( NULL );
    // Allocate space for date string
    char* date = (char*)malloc( 100 );
    // Format the time correctly [yyyy-mm-dd hh:mm:ss]
    strftime(date, 100, "[%F %T]", localtime(&t));
    return date;
}


