/* SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file gaussian_internal.h
    \brief Shared Gaussian message products for Brownian calculations and binary-trait dispersion. */
#ifndef CURUPIXA_GAUSSIAN_INTERNAL_H
#define CURUPIXA_GAUSSIAN_INTERNAL_H
#include <math.h>
/* Product of Gaussian messages, including flat and exact observations. Ratios
   avoid overflow from multiplying variances. Used by PCM, ASR and signal scores. */
static inline void
crpx_private_gaussian_product (double first_mean, double first_variance, double second_mean, double second_variance,
                               double *mean, double *variance)
{
  if (isinf (first_variance)) {
    *mean = second_mean;
    *variance = second_variance;
    return;
  }
  if (isinf (second_variance) || first_variance == 0) {
    *mean = first_mean;
    *variance = first_variance;
    return;
  }
  if (second_variance == 0) {
    *mean = second_mean;
    *variance = second_variance;
    return;
  }
  double ratio =
    first_variance <= second_variance ? first_variance / second_variance : second_variance / first_variance;
  double fraction = first_variance <= second_variance ? ratio / (1 + ratio) : 1 / (1 + ratio);
  *mean = first_mean * (1 - fraction) + second_mean * fraction;
  *variance = first_variance <= second_variance ? first_variance / (1 + ratio) : second_variance / (1 + ratio);
}
#endif
