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
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include "config.h"

#include <stdio.h>
#include <stdlib.h>     /* random number, searching, sorting, EXIT_SUCCESS */
#include <string.h>     /* string manipulation */
#include <stdarg.h>     /* Access to va_arg, va_list */
#include <stdint.h>     /* standard integer types (int32_t typedef etc.) */
#include <ctype.h>      /* char operation functions (e.g. isspace() ), case convertion */
#include <math.h>       /* standard math functions (e.g. exp() ) */
#include <float.h>      /* DBL_MAX_EXP, DBL_EPSILON constants (to avoid underflow etc) */
#include <time.h>       /* speed profiling(e.g. clock(), clock_gettime(), struct timespec ) */
#include <unistd.h>     /* sysconf(), getrandom() getentropy(), and getpid() function, used together with sys/types.h */
#include <sys/time.h>   /* random seed (e.g. gettimeofday(), struct timeval) */
#include <sys/times.h>  /* speed profiling in clock ticks (e.g. times() ) */ 
#include <sys/types.h>  /* pid_t for process ID, used together with unistd.h */
#include <sys/syscall.h>/* system calls like syscall(SYS_getrandom, buf, buflen, 0) for random noise */
#include <fcntl.h>      /* open() read() close() for /dev/urandom */
#include <sys/stat.h>   /* mkdir(); returns EEXIST from sys/types.h if dir already exist (as dir or not) */ 
#include <libgen.h>     /* standard XPG basename() - the one provided by string.h is a GNU extension, fails on macOSX */
#include <stdbool.h>    /* macros bool, true, and false */
#include <limits.h>     /* UINT_MAX etc */
#include <errno.h>      /* errno, perror(), strerror() */
#include <assert.h>    

#ifdef _OPENMP
#include <omp.h>         /* OpenMP parallel threading library when available */
#endif
#if defined HAVE_SSE || defined HAVE_AVX /* these are defined in config.h, not __SSE2__ from gcc */
//#if defined __SSE4_2__ || defined __AVX2__    /* we compile with -msse4.2 or -mno-sse (i.e. whole range is included or excluded) */
#include <immintrin.h>
#endif

/* Global constants */

#if defined(__GNUC__) && __GNUC__ >= 7
 #define CRPX_attribute_FALLTHROUGH __attribute__ ((fallthrough));
#else
 #define CRPX_attribute_FALLTHROUGH ((void)0);
#endif /* __GNUC__ >= 7 */

#define CRPX_MIN(x,y) (((x)<(y)) ? (x) : (y))
#define CRPX_MAX(x,y) (((x)>(y)) ? (x) : (y))
#define CRPX_MOD(a)   (((a)>0)   ? (a) :(-a))

#define CRPX_LOGLEVEL_FATAL   0U // enum{} would complain about unsigned
#define CRPX_LOGLEVEL_ERROR   1U
#define CRPX_LOGLEVEL_WARN    2U
#define CRPX_LOGLEVEL_INFO    3U
#define CRPX_LOGLEVEL_VERBOSE 4U
#define CRPX_LOGLEVEL_DEBUG   5U
#define crpx_logger_fatal(...)   crpx_logger_message(CRPX_LOGLEVEL_FATAL, __FILE__, __LINE__, __VA_ARGS__)
#define crpx_logger_error(...)   crpx_logger_message(CRPX_LOGLEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define crpx_logger_warning(...) crpx_logger_message(CRPX_LOGLEVEL_WARN, __FILE__, __LINE__, __VA_ARGS__)
#define crpx_logger_info(...)    crpx_logger_message(CRPX_LOGLEVEL_INFO, __FILE__, __LINE__, __VA_ARGS__)
#define crpx_logger_verbose(...) crpx_logger_message(CRPX_LOGLEVEL_VERBOSE, __FILE__, __LINE__, __VA_ARGS__)
#define crpx_logger_debug(...)   crpx_logger_message(CRPX_LOGLEVEL_DEBUG, __FILE__, __LINE__, __VA_ARGS__)

#define crpx_free(...) crpx_free_with_errmsg(__FILE__, __LINE__, __VA_ARGS__)
#define crpx_malloc(...) crpx_malloc_with_errmsg(__FILE__, __LINE__, __VA_ARGS__)
#define crpx_calloc(...) crpx_calloc_with_errmsg(__FILE__, __LINE__, __VA_ARGS__)
#define crpx_realloc(...) crpx_realloc_with_errmsg(__FILE__, __LINE__, __VA_ARGS__)
#define crpx_reallocarray(...) crpx_reallocarray_with_errmsg(__FILE__, __LINE__, __VA_ARGS__)

typedef struct {
  uint32_t id:16, loglevel_stderr:4, loglevel_file:4, error:2, sse:1, avx:1;
  uint64_t elapsed_time[2];
  FILE *logfile;
} crpx_global_struct, *crpx_global_t;


/*! \brief Memory-safe malloc() function, with default error message in case of failure. Variadic args start at cglobal */
void *crpx_malloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, size_t size);
/*! \brief Memory-safe calloc() function, with default error message in case of failure. Variadic args start at cglobal */
void * crpx_calloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, size_t nmemb, size_t size);

/*! \brief Memory-safe realloc() function, with default error message in case of failure. Variadic args start at cglobal.
 * In case of error, it does _not_ touch/free *ptr, but returns NULL; thus DO NOT "ptr = realloc (ptr)" but use a pivot/tmp pointer instead */
void *crpx_realloc_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, void *ptr, size_t size);
/*! \brief Memory-safe reallocarray() function, with default error message in case of failure. Variadic args start at cglobal. 
 * In case of error, it does _not_ touch/free *ptr, but returns NULL; thus DO NOT ptr = reallocarray (ptr)  but use a pivot/tmp instead */
void *crpx_reallocarray_with_errmsg (const char *c_file, const int c_line, crpx_global_t cglobal, void *ptr, size_t nmemb, size_t size);

void curupixa_fprintf_colour (FILE *stream, int regular, int colour, const char *message, const char *normaltext, ...);

/*! \brief Prints message to stderr,  and also to log file depending on global settings.
 * save to file not really thread safe, although we can have several FILE*s pointed to same file in append mode
 * prints to sdterr in colours and to log file in plain text, without compression. 
 * Called through macros as in crpx_logger_error(cglobal, "message") etc. variadic args start at cglobal */
void crpx_logger_message (uint8_t level, const char *c_file, const int c_line, crpx_global_t cglobal, const char *fmt, ...);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
