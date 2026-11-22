/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file clock.h
    \brief Convert dated nodes and branch rates to lengths without changing a tree. */
#ifndef CURUPIXA_CLOCK_H
#define CURUPIXA_CLOCK_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Ages increase into the past. rates=NULL uses mean_rate on every edge; otherwise
   rates contains absolute positive rates indexed by node. Root entries are ignored. Input and output arrays must not overlap.
   Zero durations are permitted. No tree lengths are changed. Output has count entries;
   count must equal node count, root length is zero. Output is unchanged on error. */
crpx_status crpx_clock_lengths (const crpx_tree *tree, const double *ages, const double *rates, double mean_rate,
                                double *lengths, size_t count);
/* Recover ages from time lengths and supplied leaf ages (NULL means contemporary
   tips at zero). Internal entries of tip_ages are ignored. Each parent's implied
   ages must agree within tolerance*max(1,|age|). workspace retains traversal memory.
   Output may be partly written on failure. Root stem is ignored. */
crpx_status crpx_tree_ages (const crpx_tree *tree, const double *tip_ages, double tolerance,
                            crpx_tree_workspace *workspace, double *ages, size_t count);
#ifdef __cplusplus
}
#endif
#endif
