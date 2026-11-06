/* Based on biomcmc-lib and N. Devillard's public-domain selection routines.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_QUICKSELECT_QUANTILE_H
#define CURUPIXA_QUICKSELECT_QUANTILE_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Selection rearranges values in place. NaNs are rejected; infinities are allowed. */
crpx_status crpx_select_double(double *values, size_t size, size_t rank, double *out);
/* Quantiles preserve input and use the biomcmc order statistic: min(floor(size*p), size-1).
 * Probabilities must lie in [0,1]. No interpolation is performed. */
crpx_status crpx_quantile_double(const double *values, size_t size, double probability, double *out);
crpx_status crpx_quantile_vector_double(const double *values, size_t size, const double *probabilities,
                                      size_t count, double *out);

#ifdef __cplusplus
}
#endif
#endif
