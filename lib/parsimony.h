/* Reimplementation of biomcmc-lib parsimony (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_PARSIMONY_H
#define CURUPIXA_PARSIMONY_H
#include "tree.h"
#include "alignment.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_parsimony_workspace crpx_parsimony_workspace;
crpx_status crpx_parsimony_workspace_new(crpx_parsimony_workspace **out);
void crpx_parsimony_workspace_free(crpx_parsimony_workspace *workspace);
/* Unordered, equal-cost Sankoff parsimony for 1..64 states, including hard
 * polytomies. Each leaf appears exactly once in tip_nodes. masks is pattern-major
 * [patterns*tip_count]; each nonzero mask is the set of allowed tip states.
 * weights may be NULL (all one). per_pattern may be NULL or patterns elements;
 * returns unweighted pattern scores and checked weighted total only on success.
 * Use one workspace per concurrent analysis. The workspace releases its temporary tree
 * binding before returning, so callers may free the tree while retaining the workspace. */
crpx_status crpx_parsimony_score(const crpx_tree *tree, const crpx_node_id *tip_nodes, size_t tip_count,
                                unsigned int states, size_t patterns, const uint64_t *masks,
                                const size_t *weights, crpx_parsimony_workspace *workspace,
                                uint64_t *total, uint64_t *per_pattern);
crpx_status crpx_parsimony_score_alignment(const crpx_tree *tree, const crpx_alignment *alignment,
                                          crpx_parsimony_workspace *workspace, uint64_t *score);
#ifdef __cplusplus
}
#endif
#endif
