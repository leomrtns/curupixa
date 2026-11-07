/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_CHARACTER_ANALYSIS_H
#define CURUPIXA_CHARACTER_ANALYSIS_H
#include "character_data.h"
#include "likelihood.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_character_analysis crpx_character_analysis;
/* A worker-owned reusable workspace gives likelihood and parsimony the same
 * tree/data mapping, pattern weights and error convention. No tree, data or model
 * ownership is transferred; the objects need remain alive only during a call.
 * Separate workers may share immutable data/models and a read-only tree. */
crpx_status crpx_character_analysis_new(crpx_character_analysis **out);
void crpx_character_analysis_free(crpx_character_analysis *analysis);
crpx_status crpx_character_analysis_likelihood(const crpx_tree *tree, const crpx_character_data *data,
                                              const crpx_substitution_model *model,
                                              const crpx_likelihood_options *options,
                                              crpx_character_analysis *analysis, double *log_likelihood);
/* Hard-polytomy Sankoff parsimony for any state count. NULL costs with count=0
 * means unordered unit change cost, evaluated in O(nodes*states*patterns).
 * Otherwise costs[parent_state*states+child_state] has states² finite nonnegative
 * entries and zero diagonal; asymmetric costs are permitted and rooting matters.
 * This general case costs O(edges*states²*patterns). Positive tip likelihood
 * entries define allowed states; their relative magnitudes do not change costs.
 * Returns checked weighted total only on success. Unit scores are represented
 * exactly up to 2^53; larger totals return OVERFLOW. A failed call clears the
 * pattern-score view. Workspace storage is O(nodes*states + patterns + taxa). */
crpx_status crpx_character_analysis_parsimony(const crpx_tree *tree, const crpx_character_data *data,
                                             const double *costs, size_t cost_count,
                                             crpx_character_analysis *analysis, double *score);
const double *crpx_character_analysis_pattern_scores(const crpx_character_analysis *analysis, size_t *count);
const double *crpx_character_analysis_pattern_log_likelihoods(const crpx_character_analysis *analysis,
                                                            size_t *count);
/* Borrowed workspace for likelihood partial-vector access; never free it. */
crpx_likelihood_workspace *crpx_character_analysis_likelihood_workspace(crpx_character_analysis *analysis);

#ifdef __cplusplus
}
#endif
#endif
