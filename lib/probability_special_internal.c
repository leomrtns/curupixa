/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT AND BSD-3-Clause
   The Temme expansion and inverse-normal rational approximation are adapted from
   scipy/xsf (SciPy 2024; Stephen L. Moshier 1984, 1987, 1989), BSD-3-Clause.
   See LICENSES/SciPy-BSD-3-Clause.txt. */
#include "probability_special_internal.h"
#include "prob_distribution.h"
#include "probability_coefficients_internal.h"
#include <math.h>
#include <float.h>
static double polynomial (double x, const double *c, unsigned n);
static double denominator (double x, const double *c, unsigned n);
static double temme (double a, double x, int upper);
static double complement (double log_p);

double
crpx_private_normal_logcdf (double x)
{ /* libm handles the central region; Mills' expansion retains log tails after erfc underflows. */
  if (x > 0) return log1p (-.5 * erfc (x * .7071067811865475244));
  if (x > -38) return log (.5 * erfc (-x * .7071067811865475244));
  if (isinf (x)) return -INFINITY;
  double inverse = 1 / x / x, term = 1, sum = 1;
  for (unsigned i = 1; i < 100; i++) {
    double next = term * (-(2 * (double) i - 1) * inverse);
    if (fabs (next) > fabs (term)) break;
    sum += next;
    term = next;
    if (fabs (term) < DBL_EPSILON * fabs (sum)) break;
  }
  return -.5 * x * x - log (-x) - .91893853320467274178 + log (sum);
}
double
crpx_private_normal_probability_quantile (double p)
{
  if (p > .13533528323661269189 && p < .86466471676338730811) {
    double y = p - .5, z = y * y;
    return (y + y * z * polynomial (z, ndtri_P0, 4) / denominator (z, ndtri_Q0, 8)) * 2.5066282746310005024;
  }
  return p > .5 ? -crpx_private_normal_quantile (log1p (-p)) : crpx_private_normal_quantile (log (p));
}
double
crpx_private_normal_quantile (double log_p)
{ /* SciPy/Cephes rational approximation, evaluated directly from log probability in the tails. */
  if (log_p == -INFINITY) return -INFINITY;
  if (log_p == 0) return INFINITY;
  bool reflect = log_p > -.6931471805599453;
  double tail = reflect ? log (-expm1 (log_p)) : log_p;
  double x;
  if (tail > -2) {
    double y = exp (log_p) - .5, z = y * y;
    x = (y + y * z * polynomial (z, ndtri_P0, 4) / denominator (z, ndtri_Q0, 8)) * 2.5066282746310005024;
  } else {
    double t = tail < -DBL_MAX / 2 ? sqrt (-tail) * 1.4142135623730950488 : sqrt (-2 * tail), z = 1 / t;
    x = -(t - log (t) / t -
          z * (t < 8 ? polynomial (z, ndtri_P1, 8) / denominator (z, ndtri_Q1, 8)
                     : polynomial (z, ndtri_P2, 8) / denominator (z, ndtri_Q2, 8)));
    if (t > 64 && t < 262144)
      for (unsigned i = 0; i < 2; i++)
        x -= (crpx_private_normal_logcdf (x) - tail) / (-x + 1 / (-x));
    if (reflect) x = -x;
  }
  return x;
}
double
crpx_private_gamma_factor (double a, double x)
{ /* Log(x^a exp(-x)/Gamma(a)); the deviance form avoids cancellation near the mode. */
  if (x == INFINITY) return -INFINITY;
  if (a >= 8 && fabs (x - a) < .4 * a) {
    double r = 1 / a, r2 = r * r;
    double correction =
      r * (1.0 / 12 +
           r2 * (-1.0 / 360 +
                 r2 * (1.0 / 1260 + r2 * (-1.0 / 1680 + r2 * (1.0 / 1188 + r2 * (-691.0 / 360360 + r2 / 156))))));
    return a * crpx_log1pmx ((x - a) / a) + .5 * log (a) - .91893853320467274178 - correction;
  }
  return (double) ((long double) a * log (x) - x - crpx_lgammafn (a, NULL));
}
double
crpx_private_gamma_tail (double a, double x, int upper)
{ /* Direct smaller-tail evaluation protects both ordinary and log probabilities. */
  if (x <= 0) return upper ? 0 : -INFINITY;
  if (isinf (x)) return upper ? -INFINITY : 0;
  if (a == 1) return upper ? -x : log (-expm1 (-x));
  if (a == .5)
    return upper ? .69314718055994530942 + crpx_private_normal_logcdf (-1.4142135623730950488 * sqrt (x))
                 : log (erf (sqrt (x)));
  double relative = fabs (x - a) / a;
  if ((a > 20 && a < 200 && relative < .3) || (a >= 200 && relative < .3)) return temme (a, x, upper);
  if (a < 1 && x < 1.1 && a * log (x) > -.5) {
    long double g = crpx_private_log_gamma_one_plus (a);
    long double power = 1, sum = 0;
    for (unsigned n = 1; n < 2000; n++) {
      power *= -x / n;
      long double term = power / (a + n);
      sum += term;
      if (fabsl (term) < LDBL_EPSILON * fabsl (sum)) break;
    }
    long double q =
      -expm1l ((long double) a * log (x) - g) - expl ((long double) a * log (x) - crpx_lgammafn (a, NULL)) * sum;
    double lq = (double) logl (q);
    return upper ? lq : complement (lq);
  }
  double factor = crpx_private_gamma_factor (a, x);
  if (x < a + 1) {
    double term = 1, sum = 1;
    for (unsigned n = 1; n < 100000; n++) {
      term *= x / (a + n);
      sum += term;
      if (term < DBL_EPSILON * sum) {
        double p = fmin (0, factor + log (sum) - log (a));
        return upper ? complement (p) : p;
      }
    }
  } else {
    double b = x + 1 - a, c = 1 / DBL_MIN, d = 1 / b, h = d;
    for (unsigned n = 1; n < 100000; n++) {
      double an = -(double) n * (n - a);
      b += 2;
      d = an * d + b;
      if (fabs (d) < DBL_MIN) d = DBL_MIN;
      c = b + an / c;
      if (fabs (c) < DBL_MIN) c = DBL_MIN;
      d = 1 / d;
      double change = d * c;
      h *= change;
      if (fabs (change - 1) < 4 * DBL_EPSILON) {
        double q = fmin (0, factor + log (h));
        return upper ? q : complement (q);
      }
    }
  }
  return NAN;
}

double
crpx_qnorm (double p, double mu, double sigma, bool log_p)
{
  if (isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) || isnan (mu) || isnan (sigma) || sigma < 0) return NAN;
  if (sigma == 0) return mu;
  return mu + sigma * (log_p ? crpx_private_normal_quantile (p) : crpx_private_normal_probability_quantile (p));
}

double
crpx_qgamma (double p, double alpha, double beta, bool log_p)
{ /* Solve in log length: the derivative stays finite even for subnormal quantiles. */
  double target = (isnan (p) || (log_p ? p > 0 : p < 0 || p > 1) ? NAN : log_p ? p : log (p));
  if (isnan (target) || !isfinite (alpha) || !isfinite (beta) || alpha < 0 || beta <= 0) return NAN;
  if (alpha == 0 || target == -INFINITY) return 0;
  if (alpha == 1) return -(log_p ? log (-expm1 (p)) : log1p (-p)) / beta;
  if (target == 0) return INFINITY;
  double lower_p = target;
  int upper = target > -.6931471805599453;
  if (upper) target = log (-expm1 (target));
  double x;
  if (alpha > 1) {
    double z = crpx_private_normal_quantile (lower_p), t = 1 - 1 / (9 * alpha) + z / (3 * sqrt (alpha));
    x = alpha * t * t * t;
  } else
    x = exp ((lower_p + crpx_private_log_gamma_one_plus (alpha)) / alpha);
  if (lower_p < -3 * alpha) x = exp ((lower_p + crpx_private_log_gamma_one_plus (alpha)) / alpha);
  if (!isfinite (x) || x <= 0) x = fmax (alpha, 1);
  if (crpx_private_gamma_tail (alpha, DBL_TRUE_MIN, 0) > lower_p)
    return exp ((lower_p + crpx_private_log_gamma_one_plus (alpha)) / alpha - log (beta));
  double low = 0, high = fmax (x, 1);
  for (unsigned i = 0; i < 1024; i++) {
    double f = crpx_private_gamma_tail (alpha, high, upper) - target;
    if (upper ? f <= 0 : f >= 0) break;
    low = high;
    high *= 2;
    if (!isfinite (high)) return INFINITY;
  }
  if (x < low || x > high) x = low + (high - low) / 2;
  for (unsigned i = 0; i < 160; i++) {
    double y = crpx_private_gamma_tail (alpha, x, upper), f = y - target;
    if (fabs (f) < 4e-15 * fmax (1, fabs (target))) return x / beta;
    if (upper ? f > 0 : f < 0)
      low = x;
    else
      high = x;
    if (nextafter (low, INFINITY) >= high) {
      double first = fabs (crpx_private_gamma_tail (alpha, low, upper) - target);
      double second = fabs (crpx_private_gamma_tail (alpha, high, upper) - target);
      return (first < second ? low : high) / beta;
    }
    double derivative = exp (crpx_private_gamma_factor (alpha, x) - y) * (upper ? -1 : 1);
    double step = f / derivative;
    double correction = 1 - .5 * step * (alpha - x - derivative);
    double next = x * exp (-(fabs (correction - 1) < .5 ? step / correction : step));
    if (next == x) return x / beta;
    if (!isfinite (next) || next <= low || next >= high) next = low == 0 ? x / 2 : low + (high - low) / 2;
    if (next == x) return x / beta;
    x = next;
  }
  return NAN;
}

static double
temme (double a, double x, int upper)
{ /* DLMF 8.12.3/4, SciPy's maintained 25-by-25 coefficient expansion. */
  double sigma = (x - a) / a, eta = copysign (sqrt (-2 * crpx_log1pmx (sigma)), sigma);
  int sign = upper ? 1 : -1;
  double z = sign * eta * sqrt (a), value = crpx_private_normal_logcdf (-z), powers[25] = {1}, sum = 0, inverse = 1,
         old = INFINITY;
  for (unsigned n = 1; n < 25; n++)
    powers[n] = powers[n - 1] * eta;
  for (unsigned k = 0; k < 25; k++) {
    double coefficient = igam_asymp_coeff_d[k][0];
    for (unsigned n = 1; n < 25; n++) {
      double term = igam_asymp_coeff_d[k][n] * powers[n];
      coefficient += term;
      if (fabs (term) < DBL_EPSILON * fabs (coefficient)) break;
    }
    double term = coefficient * inverse;
    if (fabs (term) > old) break;
    sum += term;
    if (fabs (term) < DBL_EPSILON * fabs (sum)) break;
    old = fabs (term);
    inverse /= a;
  }
  return value + log1p (sign * exp (-a * eta * eta / 2 - .91893853320467274178 - .5 * log (a) - value) * sum);
}
static double
complement (double p)
{
  return p > -.6931471805599453 ? log (-expm1 (p)) : log1p (-exp (p));
}
static double
polynomial (double x, const double *c, unsigned n)
{
  double value = c[0];
  for (unsigned i = 1; i <= n; i++)
    value = value * x + c[i];
  return value;
}
static double
denominator (double x, const double *c, unsigned n)
{
  double value = x + c[0];
  for (unsigned i = 1; i < n; i++)
    value = value * x + c[i];
  return value;
}

double
crpx_private_log_gamma_one_plus (double a)
{
  long double g;
  if (a < .125) {
    /* log Gamma(1+a), DLMF 5.7.3; zeta coefficients below are generated constants. */
    static const double zeta[] = {1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
                                  1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268,
                                  1.0040773561979443394, 1.0020083928260822144, 1.0009945751278180853,
                                  1.0004941886041194646, 1.0002460865533080483};
    long double power = a;
    g = -.57721566490153286061L * a;
    for (unsigned n = 2; n <= 32; n++) {
      power *= a;
      g += (n % 2 ? -1 : 1) * (n <= 12 ? zeta[n - 2] : 1.0L) * power / n;
    }
  } else
    g = crpx_lgammafn (1 + a, NULL);
  return (double) g;
}
