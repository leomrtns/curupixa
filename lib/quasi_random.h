/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: MIT
 */

/*! \file quasi_random.h
    \brief Generate low-discrepancy sequences for numerical exploration.
 */
#ifndef _curupixa_quasi_random_h_
#define _curupixa_quasi_random_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "global/global_variable.h"
#include "quasi_sequence.h"

/* Legacy context-based interface. Use quasi_sequence.h for deterministic reset, seek,
   checked dimensions and Sobol points. Legacy additive slopes are initialised randomly. Dimensions are checked at 1..2048;
   exhaustion fills the result with NAN. */

typedef struct {
  size_t size, leap, rem, iteration; /*!< leap and remainder are used by halton */
  double *r, *ko; /*!< r = vector with quasi-random, ko = korobov aux vector with weights (random generators) */
  crpx_global_t cglob;
  crpx_quasi_sequence *halton, *reverse;
} crpx_quasi_random_struct, *crpx_quasi_random_t;

crpx_quasi_random_t new_crpx_quasi_random (crpx_global_t cglob, size_t size);
void del_crpx_quasi_random (crpx_quasi_random_t q);
void crpx_quasi_random_reset (crpx_quasi_random_t q);
void crpx_quasi_random_next_korobov (crpx_quasi_random_t q);
void crpx_quasi_random_next_halton (crpx_quasi_random_t q);
void crpx_quasi_random_next_halton_original (crpx_quasi_random_t q);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
