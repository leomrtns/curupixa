/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/* Public distribution interfaces; numerical kernels are separate from sampling state. */
#include "prob_distribution.h"
#include "probability_special_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
static double log_probability (double p, bool log_p);
static double result (double value, bool log_p);
double
crpx_dexp_dt (double d, double lambda, double m, bool log_p)
{
  double prob_sum;
  if ((isnan (d)) || (isnan (lambda)) || (isnan (m)) || (d < 0.) || (!isfinite (lambda)) || (!isfinite (m)) ||
      (lambda <= 0.) || (m < 0.) || (m != floor (m)))
    return NAN;
  if (d > m || d < 0 || d != floor (d)) return log_p ? (-INFINITY) : 0.;
  prob_sum = crpx_expm1 (-lambda * (m + 1.)) / crpx_expm1 (-lambda); /* expm1(x) = (exp(x) - 1) when x ~ 0 */

  return log_p ? ((-lambda * d) - log (prob_sum)) : (exp (-lambda * d) / prob_sum);
}

/*! \brief cdf of discrete truncated exponential (d is discrete, m is maximum value): calculates P(D <= d) */
double
crpx_pexp_dt (double d, double lambda, double m, bool log_p)
{
  if ((isnan (d)) || (isnan (lambda)) || (isnan (m)) || (d < 0.) || (!isfinite (lambda)) || (!isfinite (m)) ||
      (lambda <= 0.) || (m < 0.) || (m != floor (m)))
    return NAN;
  if (d >= m) return log_p ? 0. : 1.;
  d = floor (d);

  return log_p ? (log (-crpx_expm1 (-lambda * (d + 1.))) - log (-crpx_expm1 (-lambda * (m + 1.))))
               : (crpx_expm1 (-lambda * (d + 1.)) / crpx_expm1 (-lambda * (m + 1.)));
}

/*! \brief quantile of discrete truncated exponential, that is, finds d s.t. P(D <= d) >= p */
double
crpx_qexp_dt (double p, double lambda, double m, bool log_p)
{
  double x;
  if ((isnan (p)) || (isnan (lambda)) || (isnan (m)) || (!isfinite (lambda)) || (!isfinite (m)) || (lambda <= 0.) ||
      (m < 0.) || (m != floor (m)))
    return NAN;

  if (log_p) {
    if (p > 0.) return NAN;
    if (p == 0.) return m;
    if (p == (-INFINITY)) return 0.;
  } else {
    if ((p < 0.) || (p > 1.)) return NAN;
    if (p == 1.) return m;
    if (p == 0.) return 0.;
  }
  if (log_p) p = exp (p);

  x = crpx_log1p (p * crpx_expm1 (-lambda * (m + 1.))); /* log1p(x) = log (1+x) */

  /* For tiny probabilities, subtraction can round to exactly -1. */
  return fmax (0, fmin (m, ceil (-x / lambda - 1.)));
}

double
crpx_lgammafn (double x, int *sign)
{ /* The re-entrant libm entry point never writes the process-global signgam. */
  extern double lgamma_r (double, int *);
  int local;
  double value = lgamma_r (x, &local);
  if (sign) *sign = local;
  return value;
}
double
crpx_log1pmx (double x)
{ /* Taylor terms avoid cancellation near zero. */
  if (fabs (x) > .125) return log1p (x) - x;
  double power = x * x, sum = -power / 2;
  for (unsigned n = 3; n < 40; n++) {
    power *= -x;
    double term = -power / n;
    sum += term;
    if (fabs (term) < DBL_EPSILON * fabs (sum)) break;
  }
  return sum;
}
double
crpx_dnorm (double x, double mu, double sigma, bool log_p)
{
  if (isnan (x) || isnan (mu) || isnan (sigma) || sigma < 0 || (isinf (x) && x == mu)) return NAN;
  if (sigma == 0) return x == mu ? INFINITY : result (-INFINITY, log_p);
  double z = (x - mu) / sigma;
  return log_p ? -.91893853320467274178 - .5 * z * z - log (sigma) : .39894228040143267794 * exp (-.5 * z * z) / sigma;
}
double
crpx_pnorm (double x, double mu, double sigma, bool log_p)
{
  if (isnan (x) || isnan (mu) || isnan (sigma) || sigma < 0 || (isinf (x) && x == mu)) return NAN;
  if (sigma == 0) return result (x < mu ? -INFINITY : 0, log_p);
  double z = (x - mu) / sigma;
  return log_p ? crpx_private_normal_logcdf (z) : .5 * erfc (-z * .7071067811865475244);
}

double
crpx_dlnorm (double x, double m, double s, bool log_p)
{
  if (isnan (x) || isnan (m) || isnan (s) || s <= 0) return NAN;
  if (x <= 0 || isinf (x)) return result (-INFINITY, log_p);
  return result (crpx_dnorm (log (x), m, s, true) - log (x), log_p);
}
double
crpx_plnorm (double x, double m, double s, bool log_p)
{
  if (isnan (x) || isnan (m) || isnan (s) || s <= 0) return NAN;
  return x <= 0 ? result (-INFINITY, log_p) : crpx_pnorm (log (x), m, s, log_p);
}
double
crpx_qlnorm (double p, double m, double s, bool log_p)
{
  if (isnan (m) || isnan (s) || s < 0 || isnan (log_probability (p, log_p))) return NAN;
  if (log_p ? p == -INFINITY : p == 0) return 0;
  if (log_p ? p == 0 : p == 1) return INFINITY;
  return exp (crpx_qnorm (p, m, s, log_p));
}
double
crpx_dgamma (double x, double a, double b, bool log_p)
{
  if (isnan (x) || !isfinite (a) || !isfinite (b) || a < 0 || b <= 0) return NAN;
  if (x < 0 || isinf (x)) return result (-INFINITY, log_p);
  if (a == 0) return x == 0 ? INFINITY : result (-INFINITY, log_p);
  if (x == 0) return a < 1 ? INFINITY : a == 1 ? result (log (b), log_p) : result (-INFINITY, log_p);
  if (x * b == 0) return result (a * (log (x) + log (b)) - crpx_lgammafn (a, NULL) - log (x), log_p);
  return result (crpx_private_gamma_factor (a, x * b) - log (x), log_p);
}
double
crpx_pgamma (double x, double a, double b, bool log_p)
{
  if (isnan (x) || !isfinite (a) || !isfinite (b) || a < 0 || b <= 0) return NAN;
  if (a == 0) return result (x < 0 ? -INFINITY : 0, log_p);
  if (x > 0 && x * b == 0) return result (a * (log (x) + log (b)) - crpx_private_log_gamma_one_plus (a), log_p);
  return result (crpx_private_gamma_tail (a, x * b, 0), log_p);
}

double
crpx_dpois (double x, double lambda, bool log_p)
{
  if (isnan (x) || !isfinite (lambda) || lambda < 0) return NAN;
  if (x < 0 || x != floor (x) || isinf (x)) return result (-INFINITY, log_p);
  if (lambda == 0) return result (x == 0 ? 0 : -INFINITY, log_p);
  if (x == 0) return result (-lambda, log_p);
  return result (crpx_private_gamma_factor (x, lambda) - log (x), log_p);
}
double
crpx_ppois (double x, double lambda, bool log_p)
{
  if (isnan (x) || !isfinite (lambda) || lambda < 0) return NAN;
  if (x < 0) return result (-INFINITY, log_p);
  if (isinf (x) || lambda == 0) return result (0, log_p);
  return result (crpx_private_gamma_tail (floor (x) + 1, lambda, 1), log_p);
}
double
crpx_qpois (double p, double lambda, bool log_p)
{
  double target = log_probability (p, log_p);
  if (isnan (target) || !isfinite (lambda) || lambda < 0 || lambda > 0x1p52) return NAN;
  if (target == -INFINITY || lambda == 0) return 0;
  if (target == 0) return INFINITY;
  double low = -1, high = fmax (1, ceil (lambda + sqrt (lambda) * crpx_private_normal_quantile (target)));
  while (crpx_ppois (high, lambda, true) < target)
    high = high * 2 + 1;
  while (high - low > 1) {
    double mid = floor (low + (high - low) / 2);
    if (crpx_ppois (mid, lambda, true) >= target)
      high = mid;
    else
      low = mid;
  }
  return high;
}
static double
log_probability (double p, bool log_p)
{
  return isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) ? NAN : log_p ? p : log (p);
}
static double
result (double x, bool log_p)
{
  return log_p ? x : exp (x);
}
double
crpx_gammafn (double x)
{
  return tgamma (x);
}
double
crpx_log1p (double x)
{
  return log1p (x);
}
double
crpx_expm1 (double x)
{
  return expm1 (x);
}

double
crpx_logspace_add (double log_x, double log_y)
{
  if (isnan (log_x) || isnan (log_y)) return NAN;
  if (log_x == INFINITY || log_y == INFINITY) return INFINITY;
  if (log_x == -INFINITY) return log_y;
  if (log_y == -INFINITY) return log_x;
  double maximum = fmax (log_x, log_y);
  return maximum + log1p (exp (fmin (log_x, log_y) - maximum));
}

double
crpx_logspace_sub (double log_x, double log_y)
{
  if (isnan (log_x) || isnan (log_y) || log_y > log_x || (log_x == INFINITY && log_y == INFINITY)) return NAN;
  if (log_x == log_y) return -INFINITY;
  if (log_y == -INFINITY || log_x == INFINITY) return log_x;
  double difference = log_y - log_x;
  return log_x + (difference > -0.69314718055994530942 ? log (-expm1 (difference)) : log1p (-exp (difference)));
}

crpx_status
crpx_discrete_gamma (double alpha, double beta, double *rates, size_t count)
{
  if (!rates || !count || !isfinite (alpha) || !isfinite (beta) || alpha <= 0 || beta <= 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / sizeof (double)) return CRPX_ERROR_OVERFLOW;
  if (count == 1) {
    double mean = alpha / beta;
    if (!isfinite (mean)) return CRPX_ERROR_OVERFLOW;
    rates[0] = mean;
    return CRPX_OK;
  }
  double *means = malloc (count * sizeof (double));
  if (!means) return CRPX_ERROR_ALLOCATION;
  double previous = 0;
  for (size_t i = 0; i < count; ++i) {
    double probability = 1;
    if (i + 1 < count) {
      double boundary = crpx_qgamma ((double) (i + 1) / (double) count, alpha, beta, false);
      probability = crpx_pgamma (boundary, alpha + 1, beta, false);
    }
    means[i] = (double) ((long double) (probability - previous) * count * alpha / beta);
    if (!isfinite (means[i]) || means[i] < 0) {
      free (means);
      return CRPX_ERROR_OUT_OF_RANGE;
    }
    previous = probability;
  }
  for (size_t i = 0; i < count; ++i)
    rates[i] = means[i];
  free (means);
  return CRPX_OK;
}


double
crpx_pgamma_tail (double x, double shape, double rate, bool upper, bool log_p)
{
  if (!upper) return crpx_pgamma (x, shape, rate, log_p);
  if (isnan (x) || !isfinite (shape) || !isfinite (rate) || shape < 0 || rate <= 0) return NAN;
  if (shape == 0) return result (x < 0 ? 0 : -INFINITY, log_p);
  if (x > 0 && x * rate == 0) {
    double p = shape * (log (x) + log (rate)) - crpx_private_log_gamma_one_plus (shape);
    return result (p > -.6931471805599453 ? log (-expm1 (p)) : log1p (-exp (p)), log_p);
  }
  return result (crpx_private_gamma_tail (shape, x * rate, 1), log_p);
}
