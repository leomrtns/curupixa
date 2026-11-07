/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * Reimplementation of biomcmc-lib distance-based tree construction.
 */
#ifndef CRPX_UPGMA_H
#define CRPX_UPGMA_H
#include "tree.h"
#include "distance_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Matrices are read-only, complete, non-negative and finite; label order follows matrix order.
 * NULL labels generates zero-based decimal names. Tree ownership passes to the caller.
 * UPGMA/single linkage use a nearest-neighbour chain: O(n^2) time/storage.
 * NJ/BioNJ use O(n^3) time and O(n^2) storage. Ties are deterministic.
 */
crpx_status crpx_tree_upgma(const struct crpx_distance_matrix_struct *distances,
                           const char *const *labels, crpx_tree **result);
crpx_status crpx_tree_single_linkage(const struct crpx_distance_matrix_struct *distances,
                                    const char *const *labels, crpx_tree **result);
/* preserve_negative=true retains NJ's estimated lengths. Otherwise negative limbs are set to
 * zero and their paired limb receives the pair distance where non-negative; reductions still use
 * the original estimates. Clamping can change patristic distances. The final pair distance is split
 * equally over an artificial binary root; its stem is absent.
 */
crpx_status crpx_tree_neighbour_joining(const struct crpx_distance_matrix_struct *distances,
                                       const char *const *labels, bool preserve_negative, crpx_tree **result);
crpx_status crpx_tree_bionj(const struct crpx_distance_matrix_struct *distances,
                           const char *const *labels, bool preserve_negative, crpx_tree **result);

#ifdef __cplusplus
}
#endif
#endif
