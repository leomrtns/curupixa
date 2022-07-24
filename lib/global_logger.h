/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file global_logger.h 
 *  \brief threaded-safe "global" variables (i.e. all functions carry it around); contains logger and PRNG.
 */

#ifndef _curupixa_lowlevel_h_
#define _curupixa_lowlevel_h_

#include "config.h"

#include <stdio.h>      /* [ANSI C C89] */
#include <stdlib.h>     /* random number, searching, sorting, EXIT_SUCCESS */
#include <string.h>     /* string manipulation */
#include <stdarg.h>     /* Access to va_arg, va_list */
#include <stdint.h>     /* standard integer types (int32_t typedef etc.) */
#include <ctype.h>      /* char operation functions (e.g. isspace() ), case convertion */
#include <math.h>       /* standard math functions (e.g. exp() ) */
#include <float.h>      /* DBL_MAX_EXP, DBL_EPSILON constants (to avoid underflow etc) */
#include <time.h>       /* speed profiling(e.g. clock(), clock_gettime(), struct timespec ) */
#include <unistd.h>     /* system values checking at runtime (e.g. sysconf() )  */
#include <sys/time.h>   /* random seed (e.g. gettimeofday(), struct timeval) */
#include <sys/times.h>  /* speed profiling in clock ticks (e.g. times() ) */ 
#include <sys/types.h>  /* pid_t for process ID, used together with unistd.h */
#include <unistd.h>     /* getpid() function, used together with sys/types.h */
#include <fcntl.h>      /* open() read() close() for /dev/urandom */
#include <sys/stat.h>   /* mkdir(); returns EEXIST from sys/types.h if dir already exist (as dir or not) */ 
#include <libgen.h>     /* standard XPG basename() - the one provided by string.h is a GNU extension, fails on macOSX */
#include <stdbool.h>    /* macros bool, true, and false */
#include <assert.h>    

#ifdef _OPENMP
#include <omp.h>         /* OpenMP parallel threading library when available */
#endif

/* Global constants */

#if defined(__GNUC__) && __GNUC__ >= 7
 #define attribute_FALLTHROUGH __attribute__ ((fallthrough));
#else
 #define attribute_FALLTHROUGH ((void)0);
#endif /* __GNUC__ >= 7 */

#define CRPX_MIN(x,y) (((x)<(y)) ? (x) : (y))
#define CRPX_MAX(x,y) (((x)>(y)) ? (x) : (y))
#define CRPX_MOD(a)   (((a)>0)   ? (a) :(-a))

/*! \brief Memory-safe malloc() function.
 *
 *  Allocates size bytes and returns a pointer to the allocated memory. An error message is thrown in case of failure.
 *  \param[in] size allocated size, in bytes
 *  \return pointer to newly allocated memory */
void *crpx_malloc (size_t size);

/*! \brief Memory-safe realloc() function.
 *
 * Changes the size of the memory block pointed to by ptr to size bytes. An error message is thrown in case of failure.
 * \param[in] size allocated size, in bytes
 * \param[in,out] ptr pointer to previously allocated memory
 * \return pointer to newly allocated memory */
void *crpx_realloc (void *ptr, size_t size);

/*! \brief Prints error message and quits program.
 *
 * similar to fprintf (stderr, ...), but exits after printing the message
 * \param[in] template va_list following same format as printf()
 * \result exits program */
void crpx_logger_error (const char *template, ...);

void crpx_logger_warning (const char *template, ...);
#endif

