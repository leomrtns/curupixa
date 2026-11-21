/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_TANGLEGRAM_H
#define CURUPIXA_TANGLEGRAM_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { crpx_node_id first, second; } crpx_tanglegram_link;
typedef struct crpx_tanglegram_workspace crpx_tanglegram_workspace;
typedef struct {
  uint64_t crossings;
  size_t evaluations;
  bool exact;
  const crpx_node_id *first_order, *second_order;
  size_t first_leaves, second_leaves;
} crpx_tanglegram_result;

/* Snapshot both ordered topologies and the supplied tip links. Trees may have
 * different sizes, unary nodes and multifurcations. Unlinked tips are retained.
 * Several links may share an endpoint (including gene copies mapped to one
 * species); links sharing either endpoint do not cross. Repeated identical
 * links are counted with their multiplicity. Labels and lengths are ignored.
 * Each input tree's stored root defines its allowed subtree orders. No input
 * pointers survive construction. Each concurrent search needs its own workspace. */
crpx_status crpx_tanglegram_workspace_new (const crpx_tree *first, const crpx_tree *second,
                                          const crpx_tanglegram_link *links, size_t count,
                                          crpx_tanglegram_workspace **result);
void crpx_tanglegram_workspace_free (crpx_tanglegram_workspace *workspace);
/* First improve both trees by subtree swaps, then enumerate all child orders
 * if budget remains. max_evaluations includes the initial order; zero requests
 * exhaustive search without a limit. Enumeration is factorial in node degrees
 * and exponential for binary trees. Each evaluation costs O(nodes + links*log(second_leaves)).
 * exact means exhaustive completion or the lower bound zero was attained.
 * Otherwise crossings is an upper bound on the minimum, not an exact distance.
 * Repeated calls reuse storage and restart from the snapshot. Result views are
 * valid until the next call/free. Neither topology is changed. */
crpx_status crpx_tanglegram_optimise (crpx_tanglegram_workspace *workspace,
                                     size_t max_evaluations, crpx_tanglegram_result *result);
#ifdef __cplusplus
}
#endif
#endif
