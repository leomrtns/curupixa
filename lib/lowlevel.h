/* 
 *This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
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

#ifndef _biomcmc_lowlevel_h_
#define _biomcmc_lowlevel_h_

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
//#include <sys/resource.h> // suggested by goptics (gpu), but don't seem needed
#include <libgen.h> /* standard XPG basename() - the one provided by string.h is a GNU extension, fails on macOSX*/

#ifdef _OPENMP
#include <omp.h>         /* OpenMP parallel threading library when available */
#endif

/* Global constants */

#define EXP_1 2.71828182845904523536028747135266 /* Euler's number */

#define true  1U /*!< Boolean TRUE  */
#define false 0U /*!< Boolean FALSE */

#if defined(__GNUC__) && __GNUC__ >= 7
 #define attribute_FALLTHROUGH __attribute__ ((fallthrough));
#else
 #define attribute_FALLTHROUGH ((void)0);
#endif /* __GNUC__ >= 7 */

#define BIOMCMC_MIN(x,y) (((x)<(y)) ? (x) : (y))
#define BIOMCMC_MAX(x,y) (((x)>(y)) ? (x) : (y))
#define BIOMCMC_MOD(a)   (((a)>0)   ? (a) :(-a))


/*! \brief Mnemonic for boolean (char is smaller than int) */
typedef unsigned char bool;

typedef struct hungarian_struct* hungarian;

struct hungarian_struct
{
  int **cost; /*! \brief cost matrix */
  int size,  /*! \brief assignment size. Cost is a square matrix, so size should be an overestimate where "missing" nodes are added w/ cost zero */
      initial_cost, /*! \brief sum of lowest input cost values for each column. The hungarian method rescales them so that minimum per column is zero */
      final_cost;   /*! \brief our final cost is on rescaled cost matrix, therefore to restore the "classical" optimal cost one should sum it with initial_cost */
  int *col_mate, *unchosen_row, *slack_row, *row_mate, *parent_row;  /*! \brief col_mate[row] with column match for row */
  double **dcost, initial_dcost, final_dcost; /*! \brief costs when working with float numbers instead of integers */
  double *row_dec_d, *col_inc_d, *slack_d; 
  int    *row_dec,   *col_inc,   *slack; /* aux vectors */
  bool is_double;
};

/*! \brief Memory-safe malloc() function.
 *
 *  Allocates size bytes and returns a pointer to the allocated memory. An error message is thrown in case of failure.
 *  \param[in] size allocated size, in bytes
 *  \return pointer to newly allocated memory */
void *biomcmc_malloc (size_t size);

/*! \brief Memory-safe realloc() function.
 *
 * Changes the size of the memory block pointed to by ptr to size bytes. An error message is thrown in case of failure.
 * \param[in] size allocated size, in bytes
 * \param[in,out] ptr pointer to previously allocated memory
 * \return pointer to newly allocated memory */
void *biomcmc_realloc (void *ptr, size_t size);

/*! \brief Prints error message and quits program.
 *
 * similar to fprintf (stderr, ...), but exits after printing the message
 * \param[in] template va_list following same format as printf()
 * \result exits program */
void biomcmc_error (const char *template, ...);

void biomcmc_warning (const char *template, ...);

/*! \brief Comparison between integers, doubles, etc. used by qsort() */
int compare_int_increasing (const void *a, const void *b);
int compare_int_decreasing (const void *a, const void *b);
int compare_uint64_increasing (const void *a, const void *b);
int compare_uint64_decreasing (const void *a, const void *b);
int compare_double_increasing (const void *a, const void *b);
int compare_double_decreasing (const void *a, const void *b);

/*! \brief edit distance between two sequences (slow), with option to allow one of sequences to terminate soon (o.w. global cost from end to end) */
uint32_t biomcmc_levenshtein_distance (const char *s1, uint32_t n1, const char *s2, uint32_t n2, uint32_t cost_sub, uint32_t cost_indel, bool skip_borders);

char * biomcmc_strrstr (const char *haystack, const char *needle); // find last occurence of needle
int biomcmc_length_common_prefix (const char *s1, const char *s2);

/* Hungarian method for bipartite matching (assignment) */
hungarian new_hungarian (int size, bool is_double);
void hungarian_reset (hungarian p);
void hungarian_update_cost (hungarian p, int row, int col, void *cost); /* pointer since we decide if int/double later */
void del_hungarian (hungarian p);
void hungarian_solve (hungarian p, int this_size);

#endif
