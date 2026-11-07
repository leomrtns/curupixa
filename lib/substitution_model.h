/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Analytic HKY formulas follow biomcmc-recomb and guenomu, with explicit state/ownership. */
#ifndef CURUPIXA_SUBSTITUTION_MODEL_H
#define CURUPIXA_SUBSTITUTION_MODEL_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_substitution_model crpx_substitution_model;
typedef enum { CRPX_BRANCH_FIXED, CRPX_BRANCH_SITEWISE_EXPONENTIAL } crpx_branch_mode;
/* Immutable reversible models. Lengths are expected substitutions per site:
 * stationary mean exit rate is one (the one-state model is constant).
 * Frequencies are positive finite weights, normalised on construction.
 * Poisson/equal-rates supports any positive state count; binary is states=2.
 * HKY uses A,C,G,T order, and positive transition/transversion multiplier kappa.
 * No eigensystem, matrix exponential or mutable global cache is used. */
crpx_status crpx_substitution_model_poisson(size_t states, crpx_substitution_model **out);
crpx_status crpx_substitution_model_f81(size_t states, const double *frequencies, crpx_substitution_model **out);
crpx_status crpx_substitution_model_hky(const double frequencies[4], double kappa, crpx_substitution_model **out);
void crpx_substitution_model_free(crpx_substitution_model *model);
size_t crpx_substitution_model_states(const crpx_substitution_model *model);
const double *crpx_substitution_model_frequencies(const crpx_substitution_model *model);
/* Row-major P[parent_state][child_state]. parameter is a fixed length or the
 * mean of an independent exponential prior; zero is the identity in both modes.
 * In exponential mode this averages a transition for one site, not a complete
 * alignment with one shared random edge length. Output is committed on success. */
crpx_status crpx_substitution_model_transition(const crpx_substitution_model *model, double parameter,
                                              crpx_branch_mode mode, double *matrix, size_t capacity);
/* Apply P to a nonnegative finite child vector. Input/output must not overlap.
 * Equal-rates/F81 use O(states) work and do not allocate a dense matrix.
 * HKY uses a four-state analytic matrix. */
crpx_status crpx_substitution_model_apply(const crpx_substitution_model *model, double parameter,
                                         crpx_branch_mode mode, const double *child, double *parent,
                                         size_t capacity);
/* Generator Q for inspection/testing; row-major, negative diagonal. */
crpx_status crpx_substitution_model_generator(const crpx_substitution_model *model, double *matrix, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
