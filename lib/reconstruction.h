/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_RECONSTRUCTION_H
#define CURUPIXA_RECONSTRUCTION_H
#include "likelihood.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_reconstruction_workspace crpx_reconstruction_workspace;
typedef struct crpx_likelihood_reconstruction crpx_likelihood_reconstruction;
typedef struct crpx_parsimony_reconstruction crpx_parsimony_reconstruction;
typedef enum {
  CRPX_PARSIMONY_UNORDERED,
  CRPX_PARSIMONY_ORDERED,
  CRPX_PARSIMONY_CUSTOM
} crpx_parsimony_cost_mode;
/* Optional, worker-owned reconstruction buffers. The constructor allocates no
 * node/state arrays. Calls retain capacities but no input identities/pointers.
 * Results own their arrays independently of workspace, tree, data and model.
 * Separate workers can share immutable inputs; synchronise tree mutation. */
crpx_status crpx_reconstruction_workspace_new(crpx_reconstruction_workspace **out);
/* Release all retained arrays and nested workspaces; owned results survive.
 * The workspace can be reused and lazily allocates again on the next call. */
void crpx_reconstruction_workspace_clear(crpx_reconstruction_workspace *workspace);
void crpx_reconstruction_workspace_free(crpx_reconstruction_workspace *workspace);
/* Full marginal posterior of each node conditional on ALL tip observations,
 * using the model's stationary root prior and the same branch/root conventions
 * as likelihood_evaluate(). options.retain_partials is ignored: reconstruction
 * needs partials. The score is the identical weighted log likelihood.
 * An impossible pattern has no conditional distribution and returns
 * INVALID_ARGUMENT, even if its weight is zero. Polytomies are supported.
 * These are marginal probabilities, not a joint MAP assignment or samples.
 * Poisson/F81 passes cost O(patterns*nodes*states), HKY has four states.
 * Full output and retained downward partials cost O(patterns*nodes*states);
 * additional worker scratch costs O(nodes*states). Constructors set *out=NULL
 * on failure; earlier owned results remain valid through later calls. */
crpx_status crpx_likelihood_reconstruct(const crpx_tree *tree, const crpx_character_data *data,
                                       const crpx_substitution_model *model,
                                       const crpx_likelihood_options *options,
                                       crpx_reconstruction_workspace *workspace,
                                       crpx_likelihood_reconstruction **out);
void crpx_likelihood_reconstruction_free(crpx_likelihood_reconstruction *result);
size_t crpx_likelihood_reconstruction_nodes(const crpx_likelihood_reconstruction *result);
size_t crpx_likelihood_reconstruction_patterns(const crpx_likelihood_reconstruction *result);
size_t crpx_likelihood_reconstruction_states(const crpx_likelihood_reconstruction *result);
double crpx_likelihood_reconstruction_score(const crpx_likelihood_reconstruction *result);
/* Read-only states-vector in character_data/model order, or NULL for invalid
 * indices. Valid until freeing this result. Sites map through site_patterns(). */
const double *crpx_likelihood_reconstruction_probabilities(const crpx_likelihood_reconstruction *result,
                                                          crpx_node_id node, size_t pattern);
/* Hard-polytomy Sankoff reconstruction, with one coherent optimal assignment
 * (smallest state on exact local ties) and all globally admissible node states.
 * UNORDERED is unit change cost, ORDERED is |parent_state-child_state|; both
 * run in O(patterns*nodes*states). CUSTOM accepts states^2 finite nonnegative
 * row-major directed costs, zero diagonal, and costs O(patterns*edges*states^2).
 * Non-CUSTOM requires costs=NULL,cost_count=0. Branch lengths are ignored.
 * Positive observation weights mean allowed states, as in parsimony scoring.
 * Reconstructions minimise each unweighted pattern separately; pattern weights
 * affect the reported total only, including when a pattern has weight zero.
 * relative_tolerance >=0 controls ONLY the ambiguity-set comparison:
 * |conditional_cost-optimum| <= tolerance*max(1,|conditional_cost|,|optimum|).
 * Zero gives exact floating-point equality. For fractional costs, choose a
 * tolerance reflecting accumulated round-off; integer costs need none.
 * Weighted unit/ordered scores are exact up to 2^53; larger totals fail.
 * This opt-in function leaves score-only API costs unchanged. */
crpx_status crpx_parsimony_reconstruct(const crpx_tree *tree, const crpx_character_data *data,
                                      crpx_parsimony_cost_mode mode, const double *costs, size_t cost_count,
                                      double relative_tolerance, crpx_reconstruction_workspace *workspace,
                                      crpx_parsimony_reconstruction **out);
void crpx_parsimony_reconstruction_free(crpx_parsimony_reconstruction *result);
size_t crpx_parsimony_reconstruction_nodes(const crpx_parsimony_reconstruction *result);
size_t crpx_parsimony_reconstruction_patterns(const crpx_parsimony_reconstruction *result);
size_t crpx_parsimony_reconstruction_states(const crpx_parsimony_reconstruction *result);
double crpx_parsimony_reconstruction_score(const crpx_parsimony_reconstruction *result);
const double *crpx_parsimony_reconstruction_pattern_scores(const crpx_parsimony_reconstruction *result);
/* Returns SIZE_MAX for an invalid index. Assignment is jointly optimal;
 * choosing states independently from optimal_states() need not be optimal. */
size_t crpx_parsimony_reconstruction_assignment(const crpx_parsimony_reconstruction *result,
                                               crpx_node_id node, size_t pattern);
const unsigned char *crpx_parsimony_reconstruction_optimal_states(const crpx_parsimony_reconstruction *result,
                                                                crpx_node_id node, size_t pattern);
#ifdef __cplusplus
}
#endif
#endif
