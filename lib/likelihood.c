/* Copyright (C) 2006-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Pruning follows biomcmc-recomb/guenomu, using explicit worker state and log-space scaling. */
#include "likelihood.h"
#include "substitution_model_internal.h"
#include "likelihood_kernel_internal.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct crpx_likelihood_workspace {
  double *partials, *scales, *pattern_logs, *scratch, *projected, *log_pi;
  crpx_log_kernel *kernels;
  size_t kernel_capacity, pi_capacity;
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
void crpx_likelihood_workspace_clear(crpx_likelihood_workspace *w)
{
  if (!w) return;
  free(w->partials); free(w->scales); free(w->pattern_logs); free(w->scratch); free(w->projected);
  free(w->log_pi); free(w->kernels); memset(w,0,sizeof(*w));
}
void crpx_likelihood_workspace_free(crpx_likelihood_workspace *w)
{ if (w) { crpx_likelihood_workspace_clear(w); free(w); } }
size_t crpx_likelihood_workspace_storage_size(const crpx_likelihood_workspace *w)
{
  return w ? sizeof(*w)+(w->partial_capacity+w->scale_capacity+w->pattern_capacity+2*w->state_capacity+
    w->pi_capacity)*sizeof(double)+w->kernel_capacity*sizeof(crpx_log_kernel) : 0;
}
struct crpx_likelihood_binding {
  crpx_tree *tree;
  const crpx_character_data *data;
  crpx_tree_workspace *traversal;
  size_t *tip_index;
};
void crpx_likelihood_binding_free(crpx_likelihood_binding *b)
{
  if (!b) return;
  crpx_tree_workspace_free(b->traversal); crpx_tree_free(b->tree);
  crpx_character_data_free((crpx_character_data *)b->data); free(b->tip_index); free(b);
}
crpx_status crpx_likelihood_binding_new(const crpx_tree *tree,const crpx_character_data *data,
                                       crpx_likelihood_binding **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  if (!tree || !data || !crpx_tree_node_count(tree)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t nodes=crpx_tree_node_count(tree),taxa=crpx_character_data_taxa(data);
  if (nodes>SIZE_MAX/sizeof(size_t) || taxa>SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_likelihood_binding *b=calloc(1,sizeof(*b));
  crpx_node_id *tips=malloc(taxa*sizeof(*tips));
  if (!b || !tips) { free(b); free(tips); return CRPX_ERROR_ALLOCATION; }
  crpx_status status=crpx_character_data_retain(data);
  if (status==CRPX_OK) { b->data=data; status=crpx_tree_clone(tree,&b->tree); }
  if (status==CRPX_OK) {
    b->tip_index=malloc(nodes*sizeof(*b->tip_index));
    if (!b->tip_index) status=CRPX_ERROR_ALLOCATION;
  }
  if (status==CRPX_OK) status=crpx_character_data_map_tips(data,b->tree,tips,taxa);
  if (status==CRPX_OK) {
    for (size_t i=0;i<nodes;++i) b->tip_index[i]=SIZE_MAX;
    for (size_t i=0;i<taxa;++i) b->tip_index[tips[i]]=i;
    status=crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&b->traversal);
  }
  if (status==CRPX_OK) status=crpx_tree_workspace_refresh(b->traversal,b->tree);
  free(tips);
  if (status!=CRPX_OK) crpx_likelihood_binding_free(b); else *out=b;
  return status;
}
static crpx_status reserve(double **array, size_t *capacity, size_t count)
{
  if (count <= *capacity) return CRPX_OK;
  if (count > SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  double *grown = realloc(*array,count*sizeof(double));
  if (!grown) return CRPX_ERROR_ALLOCATION;
  *array = grown; *capacity = count; return CRPX_OK;
}
static crpx_status evaluate_core(const crpx_tree *tree, const crpx_character_data *data,
                                    const crpx_substitution_model *m, const crpx_likelihood_options *options,
                                    crpx_likelihood_workspace *w, double *result, const crpx_likelihood_binding *binding)
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
  if (nodes > SIZE_MAX/sizeof(crpx_log_kernel) || nodes > SIZE_MAX/sizeof(size_t) ||
      taxa > SIZE_MAX/sizeof(crpx_node_id) || states > SIZE_MAX/nodes ||
      (blocks && nodes*states > SIZE_MAX/blocks) || (blocks && nodes > SIZE_MAX/blocks)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = CRPX_OK;
  crpx_tree_workspace *traversal = binding ? binding->traversal : NULL;
  crpx_node_id *tips = binding ? NULL : malloc(taxa*sizeof(*tips));
  size_t *tip_index = binding ? binding->tip_index : malloc(nodes*sizeof(*tip_index));
  if (!tip_index || (!binding && !tips)) { status=CRPX_ERROR_ALLOCATION; goto done; }
  if (nodes>w->kernel_capacity) {
    crpx_log_kernel *grown=realloc(w->kernels,nodes*sizeof(*grown));
    if (!grown) { status=CRPX_ERROR_ALLOCATION; goto done; }
    w->kernels=grown; w->kernel_capacity=nodes;
  }
  status=reserve(&w->log_pi,&w->pi_capacity,states);
  if (status!=CRPX_OK) goto done;
  crpx_log_kernel *kernels=w->kernels;
  double *log_pi=w->log_pi;
  if (!binding) {
    status = crpx_character_data_map_tips(data,tree,tips,taxa);
    if (status != CRPX_OK) goto done;
    for (size_t i = 0; i < nodes; ++i) tip_index[i] = SIZE_MAX;
    for (size_t i = 0; i < taxa; ++i) tip_index[tips[i]] = i;
    status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&traversal);
    if (status == CRPX_OK) status = crpx_tree_workspace_refresh(traversal,tree);
    if (status != CRPX_OK) goto done;
  }
  for (size_t i = 0; i < states; ++i) log_pi[i] = log(m->pi[i]);
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
    crpx_private_likelihood_prepare_kernel(m,value,o->branch_mode,kernels+node);
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
          crpx_private_likelihood_project(m,kernels+child,vectors+(size_t)child*states,log_pi,w->projected);
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
    for (size_t s = 0; s < states; ++s) log_probability = crpx_private_likelihood_log_add(log_probability,log_pi[s]+vectors[(size_t)root*states+s]);
    if (log_probability != -INFINITY) log_probability += scales[root];
    w->pattern_logs[p] = log_probability;
    if (weights[p]) total += (long double)weights[p]*log_probability;
  }
  if (isnan(total) || total > DBL_MAX || (isfinite(total) && total < -DBL_MAX)) { status = CRPX_ERROR_OVERFLOW; goto done; }
  w->nodes = nodes; w->states = states; w->patterns = patterns; w->retained = o->retain_partials; w->valid = true;
  *result = (double)total;
done:
  if (!binding) { free(tips); free(tip_index); crpx_tree_workspace_free(traversal); }
  return status;
}
crpx_status crpx_likelihood_evaluate(const crpx_tree *tree,const crpx_character_data *data,
                                    const crpx_substitution_model *model,const crpx_likelihood_options *options,
                                    crpx_likelihood_workspace *workspace,double *score)
{ return evaluate_core(tree,data,model,options,workspace,score,NULL); }
crpx_status crpx_likelihood_evaluate_prepared(const crpx_likelihood_binding *binding,
                                             const crpx_substitution_model *model,
                                             const crpx_likelihood_options *options,
                                             crpx_likelihood_workspace *workspace,double *score)
{
  if (!binding) { if (workspace) workspace->valid=false; return CRPX_ERROR_INVALID_ARGUMENT; }
  return evaluate_core(binding->tree,binding->data,model,options,workspace,score,binding);
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
