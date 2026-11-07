/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_LIKELIHOOD_H
#define CURUPIXA_LIKELIHOOD_H
#include <stdbool.h>
#include "character_data.h"
#include "substitution_model.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_likelihood_workspace crpx_likelihood_workspace;
typedef struct crpx_likelihood_binding crpx_likelihood_binding;
typedef struct {
  crpx_branch_mode branch_mode;
  const double *edge_values;
  size_t edge_count;
  double mean_branch_length;
  bool collapse_binary_root;
  double root_edge_value;
  bool retain_partials;
} crpx_likelihood_options;
/* Initialise fixed-length mode, no override array/root collapse, retained partials;
 * exponential mean/root value defaults to 0.1. NULL options means these defaults. */
void crpx_likelihood_options_init(crpx_likelihood_options *options);
crpx_status crpx_likelihood_workspace_new(crpx_likelihood_workspace **out);
void crpx_likelihood_workspace_free(crpx_likelihood_workspace *workspace);
/* Release retained numerical capacities and invalidate borrowed result views. */
void crpx_likelihood_workspace_clear(crpx_likelihood_workspace *workspace);
size_t crpx_likelihood_workspace_storage_size(const crpx_likelihood_workspace *workspace);
/* Explicit immutable preparation for repeated evaluation. Owns a snapshot of
 * tree connectivity, labels and lengths, and retains immutable data without
 * copying observations. Original tree/data may be changed/freed afterwards.
 * Original-tree revisions deliberately do not change this snapshot: construct
 * a new binding after changes. No pointer-identity caches or node references.
 * Safe to share one binding read-only across separate worker workspaces. */
crpx_status crpx_likelihood_binding_new(const crpx_tree *tree, const crpx_character_data *data,
                                       crpx_likelihood_binding **out);
void crpx_likelihood_binding_free(crpx_likelihood_binding *binding);
/* Same numerical/options contract as evaluate(); edge overrides use snapshot
 * node IDs. After warm-up for fixed dimensions, evaluation performs no mapping,
 * tree validation, traversal construction or capacity allocations. Changing
 * models/edge values rebuilds kernels, never stale model-pointer caches. */
crpx_status crpx_likelihood_evaluate_prepared(const crpx_likelihood_binding *binding,
                                             const crpx_substitution_model *model,
                                             const crpx_likelihood_options *options,
                                             crpx_likelihood_workspace *workspace, double *log_likelihood);
/* Standard pruning on any rooted tree, including unary nodes and polytomies.
 * Fixed mode uses tree edge lengths unless edge_values supplies one value per
 * node (root entry ignored). Missing/negative/nonfinite lengths return an error.
 * Exponential mode uses per-node prior means, or mean_branch_length for every
 * edge. It integrates lengths independently per EDGE AND SITE, as in the two
 * reference programs; it is not a joint shared-edge alignment integral.
 * collapse_binary_root requires two root children, replaces their transitions
 * by one root_edge_value kernel and one identity kernel, and ignores their two
 * original edge parameters. For reversible models this treats the artificial
 * root as one unrooted edge. Root stem lengths are always ignored.
 * All other degree-two vertices remain explicit model edges.
 * The workspace borrows no inputs between calls. Share immutable data/models
 * and read-only trees using one workspace per concurrent worker.
 * Log-space partials prevent product underflow. Zero-probability data return
 * OK with -INFINITY; zero-weight patterns do not contribute to the total.
 * On failure the scalar output is unchanged and cached views become invalid. */
crpx_status crpx_likelihood_evaluate(const crpx_tree *tree, const crpx_character_data *data,
                                    const crpx_substitution_model *model, const crpx_likelihood_options *options,
                                    crpx_likelihood_workspace *workspace, double *log_likelihood);
const double *crpx_likelihood_workspace_pattern_log_likelihoods(const crpx_likelihood_workspace *workspace,
                                                               size_t *count);
/* Returned log partials are normalised by max(log partial)=0; add log_scale to
 * recover their original logs. All -INFINITY denotes an impossible subtree.
 * Views remain valid until the next evaluation/free. Retention costs
 * O(patterns*nodes*states); retain_partials=false uses O(nodes*states) working
 * storage on a fresh workspace and makes partial queries return NOT_FOUND.
 * Allocated capacities are retained for reuse until workspace destruction. */
crpx_status crpx_likelihood_workspace_log_partial(const crpx_likelihood_workspace *workspace,
                                                 crpx_node_id node, size_t pattern,
                                                 const double **values, double *log_scale);
/* Scaled ordinary values exp(log partial), with the same log_scale. The borrowed
 * scratch vector is overwritten by the next ordinary-partial query/evaluation.
 * Use log_partial to retain extremely small state contributions. */
crpx_status crpx_likelihood_workspace_partial(crpx_likelihood_workspace *workspace, crpx_node_id node,
                                             size_t pattern, const double **values, double *log_scale);
#ifdef __cplusplus
}
#endif
#endif
