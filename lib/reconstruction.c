/* Full-data messages inspired by branchtraits/src/discrete.c, generalised to
 * curupixa model/data contracts. Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "reconstruction.h"
#include "likelihood_kernel_internal.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_reconstruction_workspace {
  crpx_likelihood_workspace *likelihood;
  crpx_tree_workspace *traversal;
  double *buffers;
  size_t *zeros;
  crpx_node_id *tips;
  crpx_log_kernel *kernels;
  size_t buffer_capacity, zero_capacity, tip_capacity, kernel_capacity;
};
struct crpx_likelihood_reconstruction { size_t nodes,patterns,states; double score,*probabilities; };
struct crpx_parsimony_reconstruction {
  size_t nodes,patterns,states;
  double score,*pattern_scores;
  size_t *assignment;
  unsigned char *optimal;
};
crpx_status crpx_reconstruction_workspace_new(crpx_reconstruction_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=calloc(1,sizeof(**out)); return *out?CRPX_OK:CRPX_ERROR_ALLOCATION;
}
void crpx_reconstruction_workspace_clear(crpx_reconstruction_workspace *w)
{
  if (!w) return;
  crpx_likelihood_workspace_free(w->likelihood); crpx_tree_workspace_free(w->traversal);
  free(w->buffers); free(w->zeros); free(w->tips); free(w->kernels); memset(w,0,sizeof(*w));
}
void crpx_reconstruction_workspace_free(crpx_reconstruction_workspace *w)
{ if (w) { crpx_reconstruction_workspace_clear(w); free(w); } }
static void *grow(void *array,size_t *capacity,size_t count,size_t unit,crpx_status *status)
{
  if (*status!=CRPX_OK || count<=*capacity) return array;
  if (count>SIZE_MAX/unit) { *status=CRPX_ERROR_OVERFLOW; return array; }
  void *p=realloc(array,count*unit);
  if (!p) { *status=CRPX_ERROR_ALLOCATION; return array; }
  *capacity=count; return p;
}
static crpx_status prepare(const crpx_tree *tree,const crpx_character_data *data,
                           crpx_reconstruction_workspace *w,bool kernels)
{
  size_t n=crpx_tree_node_count(tree),k=crpx_character_data_states(data),t=crpx_character_data_taxa(data);
  if (!n || !k || !t) return CRPX_ERROR_INVALID_ARGUMENT;
  if (k>SIZE_MAX/5 || n>SIZE_MAX/k || n*k>(SIZE_MAX-5*k)/3) return CRPX_ERROR_OVERFLOW;
  crpx_status status=CRPX_OK;
  w->buffers=grow(w->buffers,&w->buffer_capacity,3*n*k+5*k,sizeof(double),&status);
  w->zeros=grow(w->zeros,&w->zero_capacity,k,sizeof(size_t),&status);
  if (!kernels) w->tips=grow(w->tips,&w->tip_capacity,t,sizeof(crpx_node_id),&status);
  if (kernels) w->kernels=grow(w->kernels,&w->kernel_capacity,n,sizeof(crpx_log_kernel),&status);
  if (status==CRPX_OK && !w->traversal)
    status=crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&w->traversal);
  /* Likelihood evaluation already validates/maps observations. Parsimony owns
   * that setup here; avoid a redundant map in the likelihood reconstruction. */
  if (status==CRPX_OK && !kernels) status=crpx_character_data_map_tips(data,tree,w->tips,t);
  if (status==CRPX_OK) status=crpx_tree_workspace_refresh(w->traversal,tree);
  return status;
}
static void normalise(double *v,size_t k)
{
  double maximum=-INFINITY;
  for (size_t s=0;s<k;++s) if (v[s]>maximum) maximum=v[s];
  if (isfinite(maximum)) for (size_t s=0;s<k;++s) v[s]-=maximum;
}
static double log_sum(const double *v,size_t k)
{
  double value=-INFINITY;
  for (size_t s=0;s<k;++s) value=crpx_private_likelihood_log_add(value,v[s]);
  return value;
}
void crpx_likelihood_reconstruction_free(crpx_likelihood_reconstruction *r)
{ if (r) { free(r->probabilities); free(r); } }
size_t crpx_likelihood_reconstruction_nodes(const crpx_likelihood_reconstruction *r) { return r?r->nodes:0; }
size_t crpx_likelihood_reconstruction_patterns(const crpx_likelihood_reconstruction *r) { return r?r->patterns:0; }
size_t crpx_likelihood_reconstruction_states(const crpx_likelihood_reconstruction *r) { return r?r->states:0; }
double crpx_likelihood_reconstruction_score(const crpx_likelihood_reconstruction *r) { return r?r->score:NAN; }
const double *crpx_likelihood_reconstruction_probabilities(const crpx_likelihood_reconstruction *r,
                                                          crpx_node_id node,size_t pattern)
{ return r && node<r->nodes && pattern<r->patterns?r->probabilities+(pattern*r->nodes+node)*r->states:NULL; }
crpx_status crpx_likelihood_reconstruct(const crpx_tree *tree,const crpx_character_data *data,
                                       const crpx_substitution_model *m,const crpx_likelihood_options *options,
                                       crpx_reconstruction_workspace *w,crpx_likelihood_reconstruction **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  if (!tree || !data || !m || !w) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n=crpx_tree_node_count(tree),k=crpx_substitution_model_states(m),p=crpx_character_data_patterns(data);
  if (!n || !k || k!=crpx_character_data_states(data)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n>SIZE_MAX/k || (p && n*k>SIZE_MAX/p/sizeof(double))) return CRPX_ERROR_OVERFLOW;
  crpx_status status=prepare(tree,data,w,true);
  crpx_likelihood_reconstruction *r=NULL;
  if (status!=CRPX_OK) goto done;
  if (!w->likelihood) status=crpx_likelihood_workspace_new(&w->likelihood);
  if (status!=CRPX_OK) goto done;
  crpx_likelihood_options o;
  crpx_likelihood_options_init(&o); if (options) o=*options; o.retain_partials=true;
  double score;
  status=crpx_likelihood_evaluate(tree,data,m,&o,w->likelihood,&score);
  if (status!=CRPX_OK) goto done;
  r=calloc(1,sizeof(*r)); if (!r) { status=CRPX_ERROR_ALLOCATION; goto done; }
  r->nodes=n; r->states=k; r->patterns=p; r->score=score;
  r->probabilities=p?malloc(p*n*k*sizeof(double)):NULL;
  if (p && !r->probabilities) { status=CRPX_ERROR_ALLOCATION; goto done; }
  double *down=w->buffers,*up=down+n*k,*message=up+n*k,*scratch=message+n*k;
  double *total=scratch+k,*exclude=total+k,*temp=exclude+k,*log_pi=temp+k;
  const double *pi=crpx_substitution_model_frequencies(m);
  for (size_t s=0;s<k;++s) log_pi[s]=log(pi[s]);
  crpx_node_id root=crpx_tree_root(tree),first=CRPX_NODE_NONE,second=CRPX_NODE_NONE;
  if (o.collapse_binary_root) { first=crpx_tree_first_child(tree,root); second=crpx_tree_next_sibling(tree,first); }
  for (size_t node=0;node<n;++node) if (node!=root) {
    double value;
    if (node==first) value=o.root_edge_value;
    else if (node==second) value=0;
    else if (o.edge_values) value=o.edge_values[node];
    else if (o.branch_mode==CRPX_BRANCH_SITEWISE_EXPONENTIAL) value=o.mean_branch_length;
    else { status=crpx_tree_branch_length(tree,(crpx_node_id)node,&value); if (status!=CRPX_OK) goto done; }
    crpx_private_likelihood_prepare_kernel(m,value,o.branch_mode,w->kernels+node);
  }
  const crpx_node_id *pre=crpx_tree_workspace_preorder(w->traversal);
  const double *pattern_logs=crpx_likelihood_workspace_pattern_log_likelihoods(w->likelihood,NULL);
  for (size_t pattern=0;pattern<p;++pattern) {
    if (!isfinite(pattern_logs[pattern])) { status=CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    for (size_t node=0;node<n;++node) {
      const double *partial; double scale;
      status=crpx_likelihood_workspace_log_partial(w->likelihood,(crpx_node_id)node,pattern,&partial,&scale);
      if (status!=CRPX_OK) goto done;
      memcpy(down+node*k,partial,k*sizeof(double));
      if (node!=root) {
        crpx_private_likelihood_project(m,w->kernels+node,partial,log_pi,message+node*k);
        normalise(message+node*k,k);
      }
    }
    memcpy(up+(size_t)root*k,log_pi,k*sizeof(double));
    for (size_t pos=0;pos<n;++pos) {
      crpx_node_id node=pre[pos];
      for (size_t s=0;s<k;++s) { total[s]=0; w->zeros[s]=0; }
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child))
        for (size_t s=0;s<k;++s) {
          double value=message[(size_t)child*k+s];
          if (value==-INFINITY) ++w->zeros[s]; else total[s]+=value;
        }
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child)) {
        for (size_t s=0;s<k;++s) {
          double own=message[(size_t)child*k+s];
          size_t zeros=w->zeros[s]-(own==-INFINITY);
          exclude[s]=zeros?-INFINITY:up[(size_t)node*k+s]+total[s]-(isfinite(own)?own:0);
        }
        crpx_private_likelihood_project_transpose(m,w->kernels+child,exclude,log_pi,temp,up+(size_t)child*k);
        normalise(up+(size_t)child*k,k);
      }
      for (size_t s=0;s<k;++s) scratch[s]=down[(size_t)node*k+s]+up[(size_t)node*k+s];
      double norm=log_sum(scratch,k);
      if (!isfinite(norm)) { status=CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      double *prob=r->probabilities+(pattern*n+node)*k;
      for (size_t s=0;s<k;++s) prob[s]=exp(scratch[s]-norm);
    }
  }
  *out=r; r=NULL;
done:
  crpx_likelihood_reconstruction_free(r); crpx_tree_workspace_reset(w->traversal); return status;
}
void crpx_parsimony_reconstruction_free(crpx_parsimony_reconstruction *r)
{ if (r) { free(r->assignment); free(r->optimal); free(r->pattern_scores); free(r); } }
size_t crpx_parsimony_reconstruction_nodes(const crpx_parsimony_reconstruction *r) { return r?r->nodes:0; }
size_t crpx_parsimony_reconstruction_patterns(const crpx_parsimony_reconstruction *r) { return r?r->patterns:0; }
size_t crpx_parsimony_reconstruction_states(const crpx_parsimony_reconstruction *r) { return r?r->states:0; }
double crpx_parsimony_reconstruction_score(const crpx_parsimony_reconstruction *r) { return r?r->score:NAN; }
const double *crpx_parsimony_reconstruction_pattern_scores(const crpx_parsimony_reconstruction *r)
{ return r?r->pattern_scores:NULL; }
size_t crpx_parsimony_reconstruction_assignment(const crpx_parsimony_reconstruction *r,crpx_node_id node,size_t p)
{ return r && node<r->nodes && p<r->patterns?r->assignment[p*r->nodes+node]:SIZE_MAX; }
const unsigned char *crpx_parsimony_reconstruction_optimal_states(const crpx_parsimony_reconstruction *r,
                                                                crpx_node_id node,size_t pattern)
{ return r && node<r->nodes && pattern<r->patterns?r->optimal+(pattern*r->nodes+node)*r->states:NULL; }
static double add(double a,double b)
{ long double c=(long double)a+b; return c>DBL_MAX?INFINITY:(double)c; }
static double change(size_t a,size_t b,crpx_parsimony_cost_mode mode,const double *costs,size_t k)
{
  if (mode==CRPX_PARSIMONY_CUSTOM) return costs[a*k+b];
  if (mode==CRPX_PARSIMONY_ORDERED) return (double)(a>b?a-b:b-a);
  return a!=b;
}
static void transform(const double *input,double *out,size_t k,crpx_parsimony_cost_mode mode,
                       const double *costs,bool transpose)
{
  if (mode==CRPX_PARSIMONY_UNORDERED) {
    double minimum=INFINITY;
    for (size_t s=0;s<k;++s) if (input[s]<minimum) minimum=input[s];
    for (size_t s=0;s<k;++s) out[s]=fmin(input[s],add(minimum,1));
  } else if (mode==CRPX_PARSIMONY_ORDERED) {
    out[0]=input[0];
    for (size_t s=1;s<k;++s) out[s]=fmin(input[s],add(out[s-1],1));
    for (size_t s=k-1;s>0;--s) out[s-1]=fmin(out[s-1],add(out[s],1));
  } else for (size_t s=0;s<k;++s) {
    out[s]=INFINITY;
    for (size_t j=0;j<k;++j) {
      double v=add(input[j],transpose?costs[j*k+s]:costs[s*k+j]);
      if (v<out[s]) out[s]=v;
    }
  }
}
crpx_status crpx_parsimony_reconstruct(const crpx_tree *tree,const crpx_character_data *data,
                                      crpx_parsimony_cost_mode mode,const double *costs,size_t cost_count,
                                      double tolerance,crpx_reconstruction_workspace *w,
                                      crpx_parsimony_reconstruction **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  if (!tree || !data || !w || !isfinite(tolerance) || tolerance<0 ||
      (mode!=CRPX_PARSIMONY_UNORDERED && mode!=CRPX_PARSIMONY_ORDERED && mode!=CRPX_PARSIMONY_CUSTOM))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n=crpx_tree_node_count(tree),k=crpx_character_data_states(data),p=crpx_character_data_patterns(data);
  size_t taxa=crpx_character_data_taxa(data);
  if (!n || !k) return CRPX_ERROR_INVALID_ARGUMENT;
  if (n>SIZE_MAX/k || (p && (n>SIZE_MAX/p/sizeof(size_t) || n*k>SIZE_MAX/p)) ||
      p>SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  if (mode==CRPX_PARSIMONY_CUSTOM) {
    if (k>SIZE_MAX/k) return CRPX_ERROR_OVERFLOW;
    if (!costs || cost_count!=k*k) return CRPX_ERROR_INVALID_ARGUMENT;
    for (size_t i=0;i<cost_count;++i)
      if (!isfinite(costs[i]) || costs[i]<0 || (i/k==i%k && costs[i]!=0)) return CRPX_ERROR_INVALID_ARGUMENT;
  } else if (costs || cost_count) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status=prepare(tree,data,w,false);
  crpx_parsimony_reconstruction *r=NULL;
  if (status!=CRPX_OK) goto done;
  r=calloc(1,sizeof(*r)); if (!r) { status=CRPX_ERROR_ALLOCATION; goto done; }
  r->nodes=n; r->states=k; r->patterns=p;
  r->assignment=p?malloc(p*n*sizeof(size_t)):NULL;
  r->optimal=p?malloc(p*n*k):NULL; r->pattern_scores=p?malloc(p*sizeof(double)):NULL;
  if (p && (!r->assignment || !r->optimal || !r->pattern_scores)) { status=CRPX_ERROR_ALLOCATION; goto done; }
  double *down=w->buffers,*up=down+n*k,*message=up+n*k,*scratch=message+n*k,*total=scratch+k;
  const double *values=crpx_character_data_pattern_data(data);
  const size_t *weights=crpx_character_data_pattern_weights(data);
  const crpx_node_id *post=crpx_tree_workspace_postorder(w->traversal),*pre=crpx_tree_workspace_preorder(w->traversal);
  crpx_node_id root=crpx_tree_root(tree);
  long double weighted=0;
  for (size_t pattern=0;pattern<p;++pattern) {
    memset(down,0,n*k*sizeof(double)); memset(up,0,n*k*sizeof(double));
    for (size_t t=0;t<taxa;++t) for (size_t s=0;s<k;++s)
      down[(size_t)w->tips[t]*k+s]=values[(pattern*taxa+t)*k+s]>0?0:INFINITY;
    for (size_t pos=0;pos<n;++pos) {
      crpx_node_id node=post[pos];
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child))
        for (size_t s=0;s<k;++s) down[(size_t)node*k+s]=add(down[(size_t)node*k+s],message[(size_t)child*k+s]);
      transform(down+(size_t)node*k,message+(size_t)node*k,k,mode,costs,false);
    }
    double optimum=INFINITY;
    for (size_t s=0;s<k;++s) if (down[(size_t)root*k+s]<optimum) {
      optimum=down[(size_t)root*k+s]; r->assignment[pattern*n+root]=s;
    }
    if (!isfinite(optimum)) { status=CRPX_ERROR_OVERFLOW; goto done; }
    r->pattern_scores[pattern]=optimum; weighted+=(long double)weights[pattern]*optimum;
    if (weighted>DBL_MAX || (mode!=CRPX_PARSIMONY_CUSTOM && weighted>9007199254740992.0L)) {
      status=CRPX_ERROR_OVERFLOW; goto done;
    }
    for (size_t pos=0;pos<n;++pos) {
      crpx_node_id node=pre[pos];
      for (size_t s=0;s<k;++s) {
        double candidate=add(down[(size_t)node*k+s],up[(size_t)node*k+s]);
        double scale=fmax(1,fmax(fabs(candidate),fabs(optimum)));
        r->optimal[(pattern*n+node)*k+s]=isfinite(candidate) && fabs(candidate-optimum)<=tolerance*scale;
        total[s]=0; w->zeros[s]=0;
      }
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child))
        for (size_t s=0;s<k;++s) {
          double v=message[(size_t)child*k+s];
          if (!isfinite(v)) ++w->zeros[s]; else total[s]=add(total[s],v);
        }
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child)) {
        for (size_t s=0;s<k;++s) {
          double own=message[(size_t)child*k+s];
          size_t zeros=w->zeros[s]-(!isfinite(own));
          scratch[s]=zeros?INFINITY:add(up[(size_t)node*k+s],total[s]-(isfinite(own)?own:0));
        }
        transform(scratch,up+(size_t)child*k,k,mode,costs,true);
        size_t assigned=r->assignment[pattern*n+node],selected=0; double best=INFINITY;
        for (size_t s=0;s<k;++s) {
          double candidate=add(down[(size_t)child*k+s],change(assigned,s,mode,costs,k));
          if (candidate<best) { best=candidate; selected=s; }
        }
        r->assignment[pattern*n+child]=selected;
      }
    }
  }
  r->score=(double)weighted; *out=r; r=NULL;
done:
  crpx_parsimony_reconstruction_free(r); crpx_tree_workspace_reset(w->traversal); return status;
}
