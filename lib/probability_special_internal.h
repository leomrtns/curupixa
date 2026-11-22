/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file probability_special_internal.h
    \brief Shared log-tail incomplete gamma and inverse-normal calculations. */
#ifndef CURUPIXA_PROBABILITY_SPECIAL_INTERNAL_H
#define CURUPIXA_PROBABILITY_SPECIAL_INTERNAL_H
#if defined(__GNUC__)
#define CRPX_NUMERICAL_INTERNAL __attribute__ ((visibility ("hidden")))
#else
#define CRPX_NUMERICAL_INTERNAL
#endif
CRPX_NUMERICAL_INTERNAL double crpx_private_log_gamma_one_plus (double a);
CRPX_NUMERICAL_INTERNAL double crpx_private_gamma_tail (double a, double x, int upper);
CRPX_NUMERICAL_INTERNAL double crpx_private_normal_probability_quantile (double p);
CRPX_NUMERICAL_INTERNAL double crpx_private_normal_quantile (double log_p);
CRPX_NUMERICAL_INTERNAL double crpx_private_normal_logcdf (double x);
CRPX_NUMERICAL_INTERNAL double crpx_private_gamma_factor (double a, double x);
#endif
