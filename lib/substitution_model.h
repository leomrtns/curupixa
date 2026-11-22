/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Analytic HKY formulas follow biomcmc-recomb and guenomu, with explicit state/ownership. */
/*! \file substitution_model.h
    \brief Construct finite-state substitution models used by likelihood and ancestral reconstruction.
 */
#ifndef CURUPIXA_SUBSTITUTION_MODEL_H
#define CURUPIXA_SUBSTITUTION_MODEL_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_substitution_model crpx_substitution_model;
typedef enum { CRPX_BRANCH_FIXED, CRPX_BRANCH_SITEWISE_EXPONENTIAL } crpx_branch_mode;
/* Poisson/F81/HKY are immutable reversible models. Lengths are expected substitutions per site: stationary mean exit rate is one (the one-
   state model is constant). Frequencies are positive finite weights, normalised on construction. Poisson/equal-rates
   supports any positive state count; binary is states=2. HKY uses A,C,G,T order, and positive transition/transversion
   multiplier kappa. No eigensystem, matrix exponential or mutable global cache is used. */
crpx_status crpx_substitution_model_poisson (size_t states, crpx_substitution_model **out);
crpx_status crpx_substitution_model_f81 (size_t states, const double *frequencies, crpx_substitution_model **out);
crpx_status crpx_substitution_model_hky (const double frequencies[4], double kappa, crpx_substitution_model **out);
/* General four-state generator, including nonreversible processes. Q is row-major,
   off-diagonals>=0 and rows sum to zero. Root weights are positive and normalised;
   they need not be stationary. Rates are preserved (not rescaled to mean one).
   Fixed transitions use uniformisation with scaling/squaring; exponential-edge
   integration uses (I-mean*Q)^-1. Existing Poisson/F81/HKY paths are unchanged. */
crpx_status crpx_substitution_model_rates4 (const double generator[16], const double root_weights[4],
                                            crpx_substitution_model **out);
/* Joint states 00,01,10,11; only one trait changes in an infinitesimal event.
   rates order: 00->01,00->10,01->00,01->11,10->00,10->11,11->01,11->10.
   Equal paired rates encode independence; unconstrained rates encode dependence.
   Evaluate both through ordinary likelihood/ASR; model fitting is separate. */
crpx_status crpx_substitution_model_binary_pair (const double rates[8], const double root_weights[4],
                                                 crpx_substitution_model **out);
void crpx_substitution_model_free (crpx_substitution_model *model);
size_t crpx_substitution_model_states (const crpx_substitution_model *model);
const double *crpx_substitution_model_frequencies (const crpx_substitution_model *model);
/* Row-major P[parent_state][child_state]. parameter is a fixed length or the mean of an independent exponential prior;
   zero is the identity in both modes. In exponential mode this averages a transition for one site, not a complete
   alignment with one shared random edge length. Output is committed on success. */
crpx_status crpx_substitution_model_transition (const crpx_substitution_model *model, double parameter,
                                                crpx_branch_mode mode, double *matrix, size_t capacity);
/* Apply P to a nonnegative finite child vector. Input/output must not overlap. Equal-rates/F81 use O(states) work and
   do not allocate a dense matrix. HKY uses a four-state analytic matrix. */
crpx_status crpx_substitution_model_apply (const crpx_substitution_model *model, double parameter,
                                           crpx_branch_mode mode, const double *child, double *parent, size_t capacity);
/* Generator Q for inspection/testing; row-major, negative diagonal. */
crpx_status crpx_substitution_model_generator (const crpx_substitution_model *model, double *matrix, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
