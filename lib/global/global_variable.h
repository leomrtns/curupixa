/* 
 *This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: MIT
 */

/*! \file global_variable.h
    \brief Create and release contexts used by logging, allocation helpers and context-based random generators.
 */
#ifndef _global_global_variable_h_
#define _global_global_variable_h_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "lowlevel.h"
#include "maths_and_bits.h"
#include "hash_functions.h" // includes hash_functions_generators.h
#include "random_number.h"  // includes random_number_generators.h

/* Seed zero requests OS entropy; any other seed is reproducible. NULL level selects errors only. */
crpx_global_t crpx_global_init (uint64_t seed, const char *level_string);
void crpx_global_finalise (crpx_global_t cglob);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
