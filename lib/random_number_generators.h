/* This file is part of curupixa, a low-level library for phylogenomic analysis.
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

/*! \file random_number_generators.h
 *  \brief PRNGs   */

#ifndef _curupixa_random_number_generators_h_
#define _curupixa_random_number_generators_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "hash_functions.h" 

/* random numbers (depend on a state which is updated as new numbers are generated) */
uint64_t crpx_rng_wyhash_state64 (void *vstate);
uint64_t crpx_rng_splitmix_seed64 (void *vstate);
uint64_t crpx_lehmer_seed128 (void *vstate);
uint64_t crpx_wyrand_seed64 (void *vstate);
uint64_t crpx_rng_jenkins13_seed256 (void *vstate); // 13 bits of avalanche
uint64_t crpx_rng_jenkins19_seed256 (void *vstate); // 18.4 bits of avalanche (preferred to 13 bits)

uint64_t crpx_rng_rrmixer_seed64 (void *vstate);
uint64_t crpx_rng_moremur_seed64 (void *vstate);
uint64_t crpx_rng_romu_seed256 (void *vstate); // romu_quad: 4 x uint64_t 
uint64_t crpx_rng_romu_seed192 (void *vstate); // romu_trio: 3 x uint64_t 
uint64_t crpx_rng_romu_seed128 (void *vstate); // romu_duo: 2 x uint64_t 

uint64_t crpx_xoroshiro_pv6_seed128 (void *vstate); // 128+ V 2016
uint64_t crpx_xoroshiro_pv8_seed128 (void *vstate); // 128+ V 2018
uint64_t crpx_xoroshiro_pp_seed128 (void *vstate); // 128++
uint64_t crpx_xoroshiro_pp_seed256 (void *vstate); 
uint64_t crpx_xoroshiro_star_seed256 (void *vstate); 

uint64_t crpx_xorshift_star_seed64 (void *vstate);
uint64_t crpx_xorshift_p_seed128 (void *vstate);
uint64_t crpx_pcg_seed256 (void *vstate);

/* 32 bits */ 
uint32_t crpx_rng_abyssinian_seed128 (void *vstate); // 2 x uint64_t 
uint32_t crps_rng_widynski_seed192 (void *vstate);
uint32_t crpx_rng_jenkins8_seed128 (void *vstate); // 8 bits of avalanche (careful with similarly named 64 bits)
uint32_t crpx_rng_jenkins13_seed128 (void *vstate); // 13 bits of avalanche (careful with similarly named 64 bits)

/* seed generation */
void crpx_pcg_set_seed256 (void *vstate, uint64_t seed);
void cprx_rng_abyssinian_set_seed128 (void *vstate, uint32_t seed);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
