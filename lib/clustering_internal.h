/* SPDX-License-Identifier: MIT */
/*! \file clustering_internal.h
    Clustering requires non-negative values; the distance generator already checks finiteness. */
#ifndef CRPX_CLUSTERING_INTERNAL_H
#define CRPX_CLUSTERING_INTERNAL_H
#include "distance_generator.h"

static inline crpx_status
crpx_cluster_distance (crpx_distance_generator_t generator, size_t first, size_t second, size_t component, double *value)
{
  crpx_status status = crpx_distance_generator_get (generator, first, second, component, value);
  if (status == CRPX_OK && *value < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  return status;
}
#endif
