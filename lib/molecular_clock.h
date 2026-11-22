/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file molecular_clock.h
    \brief Sample and evaluate branch-rate models independently of an inference algorithm.
    Drummond et al. (2006), doi:10.1371/journal.pbio.0040088;
    Thorne et al. (1998), doi:10.1093/oxfordjournals.molbev.a025892. */
#ifndef CURUPIXA_MOLECULAR_CLOCK_H
#define CURUPIXA_MOLECULAR_CLOCK_H
#include "clock.h"
#include "sampling.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {
  CRPX_CLOCK_STRICT,
  CRPX_CLOCK_LOGNORMAL,
  CRPX_CLOCK_EXPONENTIAL,
  CRPX_CLOCK_AUTOCORRELATED
} crpx_clock_kind;
typedef struct {
  crpx_clock_kind kind;
  double mean_rate, dispersion;
} crpx_clock_model;
/* mean_rate is the arithmetic mean, not mean log rate. For LOGNORMAL dispersion
   is log standard deviation >0. EXPONENTIAL uses mean_rate, ignoring dispersion.
   AUTOCORRELATED: log(r_child)|r_parent ~ N(log(r_parent)-v/2,v), where
   v=dispersion^2*(duration_child+duration_parent)/2. For root children the
   ancestral rate is fixed at mean_rate at the root, and v=dispersion^2*duration_child/2.
   This is a branch-midpoint conditional-rate model; branch lengths use constant
   branch rates, not an integrated geometric Brownian rate trajectory.
   STRICT needs no random callback and fills every rate with mean_rate; its density is a point mass, 0 or -INF.
   Root rate entry always equals mean_rate and has no density factor.
   Zero-time autocorrelated transitions are point masses at the parent rate.
   Ages may include dated tips. Supply a traversal workspace per concurrent call.
   Sampling may partly write rates and consume draws on failure. */
crpx_status crpx_clock_sample (const crpx_clock_model *model, const crpx_tree *tree, const double *ages,
                               crpx_uniform_function uniform, void *state, crpx_tree_workspace *workspace,
                               double *rates, size_t count);
/* Density with respect to non-root branch rates, conditional on ages and model
   parameters. No age prior, calibration density or rate/time Jacobian is included.
   Negative/nonfinite rates are invalid; rates outside a point mass give -INF.
   Scalar output stays unchanged on error. */
crpx_status crpx_clock_log_density (const crpx_clock_model *model, const crpx_tree *tree, const double *ages,
                                    const double *rates, crpx_tree_workspace *workspace, size_t count, double *value);
#ifdef __cplusplus
}
#endif
#endif
