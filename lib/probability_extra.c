/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
#include "prob_distribution.h"
#include <math.h>
/* Exponential and chi-square use the same rate/shape kernels as gamma. */
double
crpx_dexp (double x, double rate, bool log_p)
{
  return crpx_dgamma (x, 1, rate, log_p);
}
double
crpx_pexp (double x, double rate, bool log_p)
{
  if (isnan (x) || !isfinite (rate) || rate <= 0) return NAN;
  double z = log (x) + log (rate), p = x <= 0 ? -INFINITY : z < -36 ? z : log (-expm1 (-rate * x));
  return log_p ? p : exp (p);
}
double
crpx_qexp (double p, double rate, bool log_p)
{
  if (isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) || !isfinite (rate) || rate <= 0) return NAN;
  return log_p && p < -36 ? exp (p - log (rate)) : -(log_p ? log (-expm1 (p)) : log1p (-p)) / rate;
}
double
crpx_dchisq (double x, double df, bool log_p)
{
  return df > 0 ? crpx_dgamma (x, df / 2, .5, log_p) : NAN;
}
double
crpx_pchisq (double x, double df, bool log_p)
{
  return df > 0 ? crpx_pgamma (x, df / 2, .5, log_p) : NAN;
}
double
crpx_qchisq (double p, double df, bool log_p)
{
  return df > 0 ? crpx_qgamma (p, df / 2, .5, log_p) : NAN;
}
/* Weibull uses shape and scale; its exponential special case is shape=1. */
double
crpx_dweibull (double x, double shape, double scale, bool log_p)
{
  if (isnan (x) || !isfinite (shape) || shape <= 0 || !isfinite (scale) || scale <= 0) return NAN;
  double value;
  if (x < 0 || isinf (x))
    value = -INFINITY;
  else if (x == 0)
    value = shape < 1 ? INFINITY : shape == 1 ? -log (scale) : -INFINITY;
  else {
    double z = log (x) - log (scale), power = exp (shape * z);
    value = isinf (power) ? -INFINITY : log (shape) - log (scale) + (shape - 1) * z - power;
  }
  return log_p ? value : exp (value);
}
double
crpx_pweibull (double x, double shape, double scale, bool log_p)
{
  if (isnan (x) || !isfinite (shape) || shape <= 0 || !isfinite (scale) || scale <= 0) return NAN;
  if (x <= 0) return log_p ? -INFINITY : 0;
  double z = shape * (log (x) - log (scale)), value = z < -36 ? z : log (-expm1 (-exp (z)));
  return log_p ? value : exp (value);
}
double
crpx_qweibull (double p, double shape, double scale, bool log_p)
{
  if (isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) || !isfinite (shape) || shape <= 0 || !isfinite (scale) ||
      scale <= 0)
    return NAN;
  double z = log_p && p < -36 ? p : log (-(log_p ? log (-expm1 (p)) : log1p (-p)));
  return exp (log (scale) + z / shape);
}
/* Logistic calculations use log1p/exp branches to avoid overflowing exp(-z). */
double
crpx_dlogis (double x, double location, double scale, bool log_p)
{
  if (isnan (x) || !isfinite (location) || !isfinite (scale) || scale <= 0) return NAN;
  double z = fabs ((x - location) / scale), value = -z - 2 * log1p (exp (-z)) - log (scale);
  return log_p ? value : exp (value);
}
double
crpx_plogis (double x, double location, double scale, bool log_p)
{
  if (isnan (x) || !isfinite (location) || !isfinite (scale) || scale <= 0) return NAN;
  double z = (x - location) / scale, value = z < 0 ? z - log1p (exp (z)) : -log1p (exp (-z));
  return log_p ? value : exp (value);
}
double
crpx_qlogis (double p, double location, double scale, bool log_p)
{
  if (isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) || !isfinite (location) || !isfinite (scale) || scale <= 0)
    return NAN;
  return location + scale * ((log_p ? p : log (p)) - (log_p ? log (-expm1 (p)) : log1p (-p)));
}
/* Direct upper tails avoid subtracting a small probability from a rounded one. */
double
crpx_pnorm_tail (double x, double location, double scale, bool upper, bool log_p)
{
  if (!upper) return crpx_pnorm (x, location, scale, log_p);
  if (scale == 0 && !isnan (x) && !isnan (location) && !(isinf (x) && x == location))
    return log_p ? (x < location ? 0 : -INFINITY) : (x < location ? 1 : 0);
  return crpx_pnorm (-x, -location, scale, log_p);
}
