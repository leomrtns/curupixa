/* 
 *This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.h 
 *  \brief Lowest level header file. Header file for lowlevel.c
 */

#ifndef _curupixa_h_ 
#define _curupixa_h_

#include "config.h"

#include <stdio.h>      /* [ANSI C C89] */
#include <stdlib.h>     /* random number, searching, sorting, EXIT_SUCCESS [ANSI C C89] */
#include <string.h>     /* string manipulation [ANSI C C89] */
#include <stdarg.h>     /* Access to va_arg, va_list [ANSI C C89] */
#include <stdint.h>     /* standard integer types (int32_t typedef etc.) [C99]*/
#include <ctype.h>      /* char operation functions (e.g. isspace() ), case convertion [ANSI C C89] */
#include <math.h>       /* standard math functions (e.g. exp() ) [ANSI C C89] */
#include <float.h>      /* DBL_MAX_EXP, DBL_EPSILON constants (to avoid underflow etc) */
#include <time.h>       /* speed profiling(e.g. clock(), clock_gettime(), struct timespec ) [ANSI C C89] */
#include <unistd.h>     /* system values checking at runtime (e.g. sysconf() ) [POSIX C] */
#include <sys/time.h>   /* random seed (e.g. gettimeofday(), struct timeval) [POSIX C] */
#include <sys/times.h>  /* speed profiling in clock ticks (e.g. times() ) [POSIX.1 (or GNU extension?)] */ 
#include <sys/types.h>  /* pid_t for process ID, used together with unistd.h (may not be necessary) */
#include <unistd.h>     /* getpid() function, used together with sys/types.h (may not be necessary) */
#include <fcntl.h>      /* open() read() close() for /dev/urandom */
#include <assert.h>    
#include <sys/stat.h>   /* mkdir(); returns EEXIST from sys/types.h if dir already exist (as dir or not) */ 
#include <libgen.h> /* standard XPG basename() - the one provided by string.h is a GNU extension, fails on macOSX */

#ifdef _OPENMP
#include <omp.h>         /* OpenMP parallel threading library when available */
#endif
#ifdef __SSE2__ || __SSE3__ || __SSE4_1__ || __SSE4_2__ || __AVX__ || __AVX2__ || __AVX512F__ || __AVX512DQ__ || __AVX512ER__ || __AVX512PF__
#include <immintrin.h>
#endif

#define true  1U /*!< Boolean TRUE  (also avail at stdbool.h) */
#define false 0U /*!< Boolean FALSE (also avail at stdbool.h) */

#if defined(__GNUC__) && __GNUC__ >= 7
 #define attribute_FALLTHROUGH __attribute__ ((fallthrough));
#else
 #define attribute_FALLTHROUGH ((void)0);
#endif /* __GNUC__ >= 7 */

#define CRPX_MIN(x,y) (((x)<(y)) ? (x) : (y))
#define CRPX_MAX(x,y) (((x)>(y)) ? (x) : (y))
#define CRPX_MOD(a)   (((a)>0)   ? (a) :(-a))


/*! \brief Mnemonic for boolean (char is smaller than int) */
typedef uint8_t bool;

#endif

