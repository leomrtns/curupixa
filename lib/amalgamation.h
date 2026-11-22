/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file amalgamation.h
    \brief Build a conditional-clade distribution from sampled trees; evaluate and sample amalgamated gene trees.

    Szöllősi et al. (2013), Efficient Exploration of the Space of Reconciled Gene Trees, doi:10.1093/sysbio/syt054.
    Only the empirical conditional-clade component is implemented here.
 */
#ifndef CURUPIXA_AMALGAMATION_H
#define CURUPIXA_AMALGAMATION_H
#include "tree.h"
#include "sampling.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_amalgamation crpx_amalgamation;
typedef struct crpx_amalgamation_workspace crpx_amalgamation_workspace;
typedef crpx_uniform_function crpx_amalgamation_uniform;
/* Empirical conditional clade distribution (CCD): each observed binary split receives its weighted count divided by
   the count of its parent clade. All trees must be rooted, strictly binary and contain exactly the same taxa.
   tip_nodes[tree][taxon] supplies their explicit correspondence; labels in input trees are ignored. labels supplies
   unique nonempty names for generated tips. Names may identify gene copies; species correspondence remains external.
   weights=NULL gives each sampled tree weight one. Weights must be positive. The result owns all counts, splits and
   labels, and retains no input pointers. This models a distribution on amalgamated trees, not a DTL likelihood. */
crpx_status crpx_amalgamation_new (const crpx_tree *const *trees, const crpx_node_id *const *tip_nodes,
                                   size_t tree_count, size_t taxa, const char *const *labels, const uint64_t *weights,
                                   crpx_amalgamation **out);
void crpx_amalgamation_free (crpx_amalgamation *model);
size_t crpx_amalgamation_taxa (const crpx_amalgamation *model);
size_t crpx_amalgamation_clades (const crpx_amalgamation *model);
uint64_t crpx_amalgamation_samples (const crpx_amalgamation *model);
/* One mutable workspace per concurrent evaluator. Capacity is retained; the tree is not retained. Immutable models can
   be shared between workers. */
crpx_status crpx_amalgamation_workspace_new (crpx_amalgamation_workspace **out);
void crpx_amalgamation_workspace_free (crpx_amalgamation_workspace *workspace);
void crpx_amalgamation_workspace_clear (crpx_amalgamation_workspace *workspace);
/* Product of conditional split probabilities, returned in log space. Unsupported clades or splits give OK and
   -INFINITY. Branch lengths are ignored. Binary topology is required; outputs stay unchanged on error. Preparation
   uses linear clade-table searches: O(nodes * clades * ceil(taxa/64)), rather than a hashed large-data implementation.
   Masks use O(nodes * ceil(taxa/64)) worker storage. */
crpx_status crpx_amalgamation_log_probability (const crpx_amalgamation *model, const crpx_tree *tree,
                                               const crpx_node_id *tips, crpx_amalgamation_workspace *workspace,
                                               double *log_probability);
/* Assemble a new tree from conditional splits, possibly absent from the original sample. Callback returns an
   independent uniform variate in [0,1), using only caller-owned state. Output has labels and no branch lengths. No
   global RNG, worker threads or retained mutable model state. The same CCD factors are used by probability evaluation
   and sampling. log_probability may be NULL. */
crpx_status crpx_amalgamation_sample (const crpx_amalgamation *model, crpx_amalgamation_uniform uniform, void *state,
                                      crpx_tree **out, double *log_probability);
#ifdef __cplusplus
}
#endif
#endif
