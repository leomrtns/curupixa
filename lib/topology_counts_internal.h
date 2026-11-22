/* SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file topology_counts_internal.h
    \brief Shared resolved/unresolved counting conventions for triplet and quartet comparisons. */
#ifndef CURUPIXA_TOPOLOGY_COUNTS_INTERNAL_H
#define CURUPIXA_TOPOLOGY_COUNTS_INTERNAL_H
#include "quartet_distance.h"
static inline crpx_status
crpx_private_topology_counts_add (crpx_quartet_counts *counts, unsigned first, unsigned second)
{
  if (counts->total == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  counts->total++;
  if (!first && !second)
    counts->both_unresolved++;
  else if (!first)
    counts->second_resolved_only++;
  else if (!second)
    counts->first_resolved_only++;
  else if (first == second)
    counts->resolved_agreement++;
  else
    counts->resolved_conflict++;
  return CRPX_OK;
}
#endif
