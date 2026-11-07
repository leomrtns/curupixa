/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Pruning follows biomcmc-recomb/guenomu, using explicit worker state and log-space scaling. */
#include "likelihood.h"
#include "substitution_model_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct crpx_likelihood_workspace {
  double *partials, *scales, *pattern_logs, *scratch, *projected;
  size_t partial_capacity, scale_capacity, pattern_capacity, state_capacity;
  size_t nodes, patterns, states;
  bool valid, retained;
};
void crpx_likelihood_options_init(crpx_likelihood_options *o)
{
  if (o) *o = (crpx_likelihood_options){CRPX_BRANCH_FIXED,NULL,0,0.1,false,0.1,true};
}
crpx_status crpx_likelihood_workspace_new(crpx_likelihood_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = calloc(1,sizeof(**out));
  return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void crpx_likelihood_workspace_free(crpx_likelihood_workspace *w)
{
  if (w) { free(w->partials); free(w->scales); free(w->pattern_logs); free(w->scratch); free(w->projected); free(w); }
}
static crpx_status reserve(double **array, size_t *capacity, size_t count)
{
  if (count <= *capacity) return CRPX_OK;
  if (count > SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  double *grown = realloc(*array,count*sizeof(double));
  if (!grown) return CRPX_ERROR_ALLOCATION;
  *array = grown; *capacity = count; return CRPX_OK;
}
static double log_add(double a, double b)
{
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a < b) { double swap = a; a = b; b = swap; }
  return a+log1p(exp(b-a));
}
typedef struct { double loss, factor, matrix[16]; } log_kernel;
static void prepare_log_kernel(const crpx_substitution_model *m, double t, crpx_branch_mode mode, log_kernel *k)
{
  if (m->hky) {
    crpx_private_hky_log_matrix(m,t,mode,k->matrix);
  } else if (!t || !m->base_decay) { k->factor = 0; k->loss = -INFINITY; }
  else {
    double logx = log(m->base_decay)+log(t);
    if (mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) {
      double denom = log_add(0,logx);
      k->factor = -denom; k->loss = logx-denom;
    } else {
      double x = m->base_decay*t;
      k->factor = -x;
      k->loss = x ? log(-expm1(-x)) : logx;
    }
  }
}
static void project(const crpx_substitution_model *m, const log_kernel *k, const double *child,
                    const double *log_pi, double *out)
{
  if (m->hky) {
    for (size_t i = 0; i < 4; ++i) {
      double value = -INFINITY;
      for (size_t j = 0; j < 4; ++j) value = log_add(value,k->matrix[i*4+j]+child[j]);
      out[i] = value;
    }
  } else {
    double value = -INFINITY;
    for (size_t j = 0; j < m->states; ++j) value = log_add(value,log_pi[j]+child[j]);
    for (size_t i = 0; i < m->states; ++i) out[i] = log_add(k->loss+value,k->factor+child[i]);
  }
}
crpx_status crpx_likelihood_evaluate(const crpx_tree *tree, const crpx_character_data *data,
                                    const crpx_substitution_model *m, const crpx_likelihood_options *options,
                                    crpx_likelihood_workspace *w, double *result)
{
  if (w) w->valid = false;
  if (!tree || !data || !m || !w || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_likelihood_options defaults; crpx_likelihood_options_init(&defaults);
  const crpx_likelihood_options *o = options ? options : &defaults;
  size_t nodes = crpx_tree_node_count(tree), patterns = crpx_character_data_patterns(data);
  size_t states = m->states, taxa = crpx_character_data_taxa(data);
  if (!nodes || crpx_character_data_states(data) != states ||
      (o->branch_mode != CRPX_BRANCH_FIXED && o->branch_mode != CRPX_BRANCH_SITEWISE_EXPONENTIAL) ||
      (o->edge_values ? o->edge_count != nodes : o->edge_count != 0)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t blocks = o->retain_partials ? patterns : (patterns ? 1 : 0);
  if (nodes > SIZE_MAX/sizeof(log_kernel) || nodes > SIZE_MAX/sizeof(size_t) ||
      taxa > SIZE_MAX/sizeof(crpx_node_id) || states > SIZE_MAX/nodes ||
      (blocks && nodes*states > SIZE_MAX/blocks) || (blocks && nodes > SIZE_MAX/blocks)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = CRPX_OK;
  crpx_tree_workspace *traversal = NULL;
  crpx_node_id *tips = malloc(taxa*sizeof(*tips));
  size_t *tip_index = malloc(nodes*sizeof(*tip_index));
  log_kernel *kernels = calloc(nodes,sizeof(*kernels));
  double *log_pi = malloc(states*sizeof(double));
  if (!tips || !tip_index || !kernels || !log_pi) { status = CRPX_ERROR_ALLOCATION; goto done; }
  status = crpx_character_data_map_tips(data,tree,tips,taxa);
  if (status != CRPX_OK) goto done;
  for (size_t i = 0; i < nodes; ++i) tip_index[i] = SIZE_MAX;
  for (size_t i = 0; i < taxa; ++i) tip_index[tips[i]] = i;
  for (size_t i = 0; i < states; ++i) log_pi[i] = log(m->pi[i]);
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&traversal);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh(traversal,tree);
  if (status != CRPX_OK) goto done;
  crpx_node_id root = crpx_tree_root(tree), first = CRPX_NODE_NONE, second = CRPX_NODE_NONE;
  if (o->collapse_binary_root) {
    first = crpx_tree_first_child(tree,root);
    if (first != CRPX_NODE_NONE) second = crpx_tree_next_sibling(tree,first);
    if (second == CRPX_NODE_NONE || crpx_tree_next_sibling(tree,second) != CRPX_NODE_NONE ||
        !isfinite(o->root_edge_value) || o->root_edge_value < 0) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
  }
  for (size_t node = 0; node < nodes; ++node) if (node != root) {
    double value;
    if (node == first) value = o->root_edge_value;
    else if (node == second) value = 0;
    else if (o->edge_values) value = o->edge_values[node];
    else if (o->branch_mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) value = o->mean_branch_length;
    else {
      status = crpx_tree_branch_length(tree,(crpx_node_id)node,&value);
      if (status != CRPX_OK) goto done;
    }
    if (!isfinite(value) || value < 0) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    prepare_log_kernel(m,value,o->branch_mode,kernels+node);
    if (m->hky) for (size_t entry = 0; entry < 16; ++entry)
      if (isnan(kernels[node].matrix[entry]) || kernels[node].matrix[entry] == INFINITY) {
        status = CRPX_ERROR_OVERFLOW; goto done;
      }
  }
  status = reserve(&w->partials,&w->partial_capacity,blocks*nodes*states);
  if (status == CRPX_OK) status = reserve(&w->scales,&w->scale_capacity,blocks*nodes);
  if (status == CRPX_OK) status = reserve(&w->pattern_logs,&w->pattern_capacity,patterns);
  if (status == CRPX_OK && states > w->state_capacity) {
    double *scratch = malloc(states*sizeof(double)), *projected = malloc(states*sizeof(double));
    if (!scratch || !projected) { free(scratch); free(projected); status = CRPX_ERROR_ALLOCATION; }
    else { free(w->scratch); free(w->projected); w->scratch = scratch; w->projected = projected; w->state_capacity = states; }
  }
  if (status != CRPX_OK) goto done;
  const double *observations = crpx_character_data_pattern_data(data);
  const size_t *weights = crpx_character_data_pattern_weights(data);
  const crpx_node_id *order = crpx_tree_workspace_postorder(traversal);
  long double total = 0;
  for (size_t p = 0; p < patterns; ++p) {
    size_t block = o->retain_partials ? p : 0;
    double *vectors = w->partials+block*nodes*states, *scales = w->scales+block*nodes;
    for (size_t k = 0; k < nodes; ++k) {
      crpx_node_id node = order[k], child = crpx_tree_first_child(tree,node);
      double *v = vectors+(size_t)node*states;
      scales[node] = 0;
      if (child == CRPX_NODE_NONE) {
        const double *input = observations+(p*taxa+tip_index[node])*states;
        double maximum = -INFINITY;
        for (size_t s = 0; s < states; ++s) { v[s] = log(input[s]); if (v[s] > maximum) maximum = v[s]; }
        for (size_t s = 0; s < states; ++s) v[s] -= maximum;
        scales[node] = maximum;
      } else {
        for (size_t s = 0; s < states; ++s) v[s] = 0;
        for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling(tree,child)) {
          project(m,kernels+child,vectors+(size_t)child*states,log_pi,w->projected);
          double maximum = -INFINITY;
          for (size_t s = 0; s < states; ++s) { v[s] += w->projected[s]; if (v[s] > maximum) maximum = v[s]; }
          if (maximum == -INFINITY) { scales[node] = 0; continue; }
          scales[node] += scales[child]+maximum;
          if (!isfinite(scales[node])) { status = CRPX_ERROR_OVERFLOW; goto done; }
          for (size_t s = 0; s < states; ++s) v[s] -= maximum;
        }
      }
    }
    double log_probability = -INFINITY;
    for (size_t s = 0; s < states; ++s) log_probability = log_add(log_probability,log_pi[s]+vectors[(size_t)root*states+s]);
    if (log_probability != -INFINITY) log_probability += scales[root];
    w->pattern_logs[p] = log_probability;
    if (weights[p]) total += (long double)weights[p]*log_probability;
  }
  if (isnan(total) || total > DBL_MAX || (isfinite(total) && total < -DBL_MAX)) { status = CRPX_ERROR_OVERFLOW; goto done; }
  w->nodes = nodes; w->states = states; w->patterns = patterns; w->retained = o->retain_partials; w->valid = true;
  *result = (double)total;
done:
  free(tips); free(tip_index); free(kernels); free(log_pi); crpx_tree_workspace_free(traversal);
  return status;
}
const double *crpx_likelihood_workspace_pattern_log_likelihoods(const crpx_likelihood_workspace *w, size_t *count)
{
  if (count) *count = w && w->valid ? w->patterns : 0;
  return w && w->valid ? w->pattern_logs : NULL;
}
crpx_status crpx_likelihood_workspace_log_partial(const crpx_likelihood_workspace *w, crpx_node_id node,
                                                 size_t pattern, const double **values, double *scale)
{
  if (!w || !values || !scale) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!w->valid || !w->retained) return CRPX_ERROR_NOT_FOUND;
  if (node >= w->nodes || pattern >= w->patterns) return CRPX_ERROR_OUT_OF_RANGE;
  *values = w->partials+(pattern*w->nodes+node)*w->states;
  *scale = w->scales[pattern*w->nodes+node]; return CRPX_OK;
}
crpx_status crpx_likelihood_workspace_partial(crpx_likelihood_workspace *w, crpx_node_id node, size_t pattern,
                                             const double **values, double *scale)
{
  if (!values || !scale) return CRPX_ERROR_INVALID_ARGUMENT;
  const double *logs; double offset;
  crpx_status status = crpx_likelihood_workspace_log_partial(w,node,pattern,&logs,&offset);
  if (status != CRPX_OK) return status;
  for (size_t s = 0; s < w->states; ++s) w->scratch[s] = exp(logs[s]);
  *values = w->scratch; *scale = offset; return CRPX_OK;
}
