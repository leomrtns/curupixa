/* Derived from guenomu/genefam-dist topology_splitset (Leonardo de Oliveira Martins, 2009).
 * Copyright (C) 2009-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/*! \file spr_distance.h
    \brief Calculate the historical greedy dSPR approximation from ordinary or previously prepared splitsets.
 */
#ifndef CURUPIXA_SPR_DISTANCE_H
#define CURUPIXA_SPR_DISTANCE_H
#include "splitset_distances.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_spr_workspace crpx_spr_workspace;
typedef struct {
  size_t rf, spr, spr_extra;
  uint64_t hdist, reduced_hdist;
} crpx_spr_result;
/* Optional worker-owned scratch: construction allocates no split/assignment arrays. Capacity is reused until
   clear/free. No borrowed input pointers survive a call. */
crpx_status crpx_spr_workspace_new (crpx_spr_workspace **out);
void crpx_spr_workspace_free (crpx_spr_workspace *workspace);
void crpx_spr_workspace_clear (crpx_spr_workspace *workspace);
size_t crpx_spr_workspace_storage_size (const crpx_spr_workspace *workspace);
/* Historical greedy dSPR heuristic, NOT exact SPR or a guaranteed bound/metric. rf is exact unrooted split symmetric
   difference. hdist/reduced_hdist use the historical zero dummy-edge assignment penalty; crpx_splitset_compare()
   instead uses a positive unmatched penalty. spr_extra is the separate crossing-prune diagnostic, not silently added
   to spr. Polytomies/mapped multi-copy species trees follow split compatibility: exhaustion of either splitset stops
   pruning. Tie choices are deterministic for a fixed taxon bit order; heuristics can vary with taxon order and
   direction. Inputs remain unchanged, result only on success. Workspace O((a+b)*ceil(taxa/64) + max(a,b)^2), not all-
   pairs bitset storage. */
crpx_status crpx_splitset_spr_distance (const crpx_splitset *a, const crpx_splitset *b, crpx_spr_workspace *workspace,
                                        crpx_spr_result *result);
/* Same unique leaf labels, independent of their node order. */
crpx_status crpx_tree_spr_distance (const crpx_tree *a, const crpx_tree *b, crpx_spr_workspace *workspace,
                                    crpx_spr_result *result);
#ifdef __cplusplus
}
#endif
#endif
