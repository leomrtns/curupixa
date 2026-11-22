/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file particle_weights.h
    \brief Normalise importance weights, choose temperatures and resample particle indices.
    Del Moral, Doucet and Jasra (2006), doi:10.1111/j.1467-9868.2006.00553.x. */
#ifndef CURUPIXA_PARTICLE_WEIGHTS_H
#define CURUPIXA_PARTICLE_WEIGHTS_H
#include "sampling.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Log weights may include -INF, but at least one must be finite. weights may alias
   log_weights. Outputs: normalised ordinary weights, log(sum exp(log_weights)),
   and ESS=1/sum(weights^2). No allocations or hidden state. */
crpx_status crpx_particle_normalise (const double *log_weights, size_t count, double *weights, double *log_sum,
                                     double *ess);
/* Largest next beta in (beta,1] with conditional ESS >= fraction*count.
   CESS=count*(sum w exp(delta*score))^2/sum w exp(2*delta*score).
   Input weights are normalised and scores finite. fraction lies in (0,1).
   Bisection controls the temperature, not the target density. */
crpx_status crpx_particle_temperature (const double *weights, const double *scores, size_t count, double beta,
                                       double fraction, double *next_beta);
/* Input weights are normalised nonnegative probabilities. systematic uses one
   draw; stratified=true uses one per output. Returned indices refer to old
   particles; copying/releasing their states belongs to the application.
   Outputs may be partly written after a bad random draw. */
crpx_status crpx_particle_resample (const double *weights, size_t count, size_t offspring, bool stratified,
                                    crpx_uniform_function uniform, void *state, size_t *indices);
#ifdef __cplusplus
}
#endif
#endif
