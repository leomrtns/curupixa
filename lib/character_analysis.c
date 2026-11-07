/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "character_analysis.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

struct crpx_character_analysis {
  crpx_likelihood_workspace *likelihood;
  crpx_tree_workspace *traversal;
  double *node_costs, *scores;
  crpx_node_id *tips;
  size_t cost_capacity, score_capacity, tip_capacity, valid_patterns;
};

crpx_status crpx_character_analysis_new(crpx_character_analysis **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_character_analysis *analysis = calloc(1,sizeof(*analysis));
  if (!analysis) return CRPX_ERROR_ALLOCATION;
  crpx_status status = crpx_likelihood_workspace_new(&analysis->likelihood);
  if (status == CRPX_OK) status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&analysis->traversal);
  if (status != CRPX_OK) { crpx_character_analysis_free(analysis); return status; }
  *out = analysis; return CRPX_OK;
}
void crpx_character_analysis_free(crpx_character_analysis *analysis)
{
  if (!analysis) return;
  crpx_likelihood_workspace_free(analysis->likelihood); crpx_tree_workspace_free(analysis->traversal);
  free(analysis->node_costs); free(analysis->scores); free(analysis->tips); free(analysis);
}
crpx_status crpx_character_analysis_likelihood(const crpx_tree *tree, const crpx_character_data *data,
                                              const crpx_substitution_model *model,
                                              const crpx_likelihood_options *options,
                                              crpx_character_analysis *analysis, double *score)
{
  if (!analysis) return CRPX_ERROR_INVALID_ARGUMENT;
  return crpx_likelihood_evaluate(tree,data,model,options,analysis->likelihood,score);
}

static crpx_status reserve_work(crpx_character_analysis *analysis, size_t nodes, size_t states,
                                size_t patterns, size_t taxa)
{
  if (nodes > SIZE_MAX/states || nodes*states > SIZE_MAX/sizeof(double) ||
      patterns > SIZE_MAX/sizeof(double) || taxa > SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  size_t costs = nodes*states;
  if (costs > analysis->cost_capacity) {
    double *memory = realloc(analysis->node_costs,costs*sizeof(double));
    if (!memory) return CRPX_ERROR_ALLOCATION;
    analysis->node_costs = memory; analysis->cost_capacity = costs;
  }
  if (patterns > analysis->score_capacity) {
    double *memory = realloc(analysis->scores,patterns*sizeof(double));
    if (!memory) return CRPX_ERROR_ALLOCATION;
    analysis->scores = memory; analysis->score_capacity = patterns;
  }
  if (taxa > analysis->tip_capacity) {
    crpx_node_id *memory = realloc(analysis->tips,taxa*sizeof(crpx_node_id));
    if (!memory) return CRPX_ERROR_ALLOCATION;
    analysis->tips = memory; analysis->tip_capacity = taxa;
  }
  return CRPX_OK;
}

static double add_cost(double a, double b)
{
  long double result = (long double)a+b;
  return result > DBL_MAX ? INFINITY : (double)result;
}

crpx_status crpx_character_analysis_parsimony(const crpx_tree *tree, const crpx_character_data *data,
                                             const double *costs, size_t cost_count,
                                             crpx_character_analysis *analysis, double *score)
{
  if (analysis) analysis->valid_patterns = 0;
  if (!tree || !data || !analysis || !score) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t nodes = crpx_tree_node_count(tree), states = crpx_character_data_states(data);
  size_t taxa = crpx_character_data_taxa(data), patterns = crpx_character_data_patterns(data);
  if (!nodes || !states || !taxa || (!costs && cost_count)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (costs) {
    if (states > SIZE_MAX/states) return CRPX_ERROR_OVERFLOW;
    if (cost_count != states*states) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t i = 0; i < cost_count; ++i)
      if (!isfinite(costs[i]) || costs[i] < 0 || (i/states == i%states && costs[i] != 0))
        return CRPX_ERROR_INVALID_ARGUMENT;
  }
  crpx_status status = reserve_work(analysis,nodes,states,patterns,taxa);
  if (status == CRPX_OK) status = crpx_character_data_map_tips(data,tree,analysis->tips,taxa);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh(analysis->traversal,tree);
  if (status != CRPX_OK) { crpx_tree_workspace_reset(analysis->traversal); return status; }
  const crpx_node_id *order = crpx_tree_workspace_postorder(analysis->traversal);
  const double *values = crpx_character_data_pattern_data(data);
  const size_t *weights = crpx_character_data_pattern_weights(data);
  long double total = 0;
  uint64_t integer_total = 0;
  for (size_t pattern = 0; pattern < patterns; ++pattern) {
    for (size_t taxon = 0; taxon < taxa; ++taxon) {
      double *tip = analysis->node_costs+(size_t)analysis->tips[taxon]*states;
      for (size_t state = 0; state < states; ++state)
        tip[state] = values[(pattern*taxa+taxon)*states+state] > 0 ? 0 : INFINITY;
    }
    for (size_t i = 0; i < nodes; ++i) {
      crpx_node_id node = order[i], child = crpx_tree_first_child(tree,node);
      if (child == CRPX_NODE_NONE) continue;
      double *parent_costs = analysis->node_costs+(size_t)node*states;
      for (size_t state = 0; state < states; ++state) parent_costs[state] = 0;
      for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling(tree,child)) {
        const double *child_costs = analysis->node_costs+(size_t)child*states;
        if (!costs) {
          double minimum = INFINITY;
          for (size_t state = 0; state < states; ++state) if (child_costs[state] < minimum) minimum = child_costs[state];
          for (size_t state = 0; state < states; ++state) {
            double alternative = add_cost(minimum,1);
            double best = child_costs[state] < alternative ? child_costs[state] : alternative;
            parent_costs[state] = add_cost(parent_costs[state],best);
          }
        } else for (size_t state = 0; state < states; ++state) {
          double minimum = INFINITY;
          for (size_t next = 0; next < states; ++next) {
            double candidate = add_cost(child_costs[next],costs[state*states+next]);
            if (candidate < minimum) minimum = candidate;
          }
          parent_costs[state] = add_cost(parent_costs[state],minimum);
        }
      }
    }
    double minimum = INFINITY;
    const double *root = analysis->node_costs+(size_t)crpx_tree_root(tree)*states;
    for (size_t state = 0; state < states; ++state) if (root[state] < minimum) minimum = root[state];
    if (!isfinite(minimum)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    if (!costs) {
      const uint64_t limit = UINT64_C(9007199254740992);
      uint64_t unit_score = (uint64_t)minimum;
      if (unit_score && (weights[pattern] > limit/unit_score ||
          integer_total > limit-unit_score*weights[pattern])) { status = CRPX_ERROR_OVERFLOW; goto done; }
      integer_total += unit_score*weights[pattern];
      total = integer_total;
    } else {
      total += (long double)minimum*weights[pattern];
      if (total > DBL_MAX) { status = CRPX_ERROR_OVERFLOW; goto done; }
    }
    analysis->scores[pattern] = minimum;
  }
  *score = (double)total; analysis->valid_patterns = patterns;
done:
  crpx_tree_workspace_reset(analysis->traversal); return status;
}

const double *crpx_character_analysis_pattern_scores(const crpx_character_analysis *analysis, size_t *count)
{
  if (count) *count = analysis ? analysis->valid_patterns : 0;
  return analysis && analysis->valid_patterns ? analysis->scores : NULL;
}
const double *crpx_character_analysis_pattern_log_likelihoods(const crpx_character_analysis *analysis, size_t *count)
{
  if (!analysis) { if (count) *count = 0; return NULL; }
  return crpx_likelihood_workspace_pattern_log_likelihoods(analysis->likelihood,count);
}
crpx_likelihood_workspace *crpx_character_analysis_likelihood_workspace(crpx_character_analysis *analysis)
{ return analysis ? analysis->likelihood : NULL; }
