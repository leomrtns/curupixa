/* Linear Gaussian pruning and PIC adapted from branchtraits/src/continuous.c.
 * Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "continuous_traits.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

#define LOG_2PI 1.837877066409345483560659472811
typedef struct {
  size_t nodes, tips;
  crpx_node_id root;
  crpx_tree_workspace *workspace;
  double *length, *values;
} continuous_tree;
struct crpx_brownian_result {
  size_t nodes, tips;
  double root_mean, sigma2, log_likelihood;
  double *mean, *variance;
};
struct crpx_pic_result {
  size_t count;
  crpx_node_id *node;
  double *contrast, *variance;
};
static void tree_clear(continuous_tree *index)
{
  crpx_tree_workspace_free(index->workspace); free(index->length); free(index->values);
}
static crpx_status tree_prepare(const crpx_tree *tree, const crpx_node_id *tips,
                                const double *values, size_t count, continuous_tree *index)
{
  if (!tree || !tips || !values) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status=crpx_tree_validate(tree);
  if (status!=CRPX_OK) return status;
  index->nodes=crpx_tree_node_count(tree); index->tips=crpx_tree_leaf_count(tree);
  index->root=crpx_tree_root(tree);
  if (!index->nodes || index->tips<2) return CRPX_ERROR_TOPOLOGY;
  if (count!=index->tips) return CRPX_ERROR_INVALID_ARGUMENT;
  if (index->nodes>SIZE_MAX/sizeof(double)) return CRPX_ERROR_OVERFLOW;
  index->length=calloc(index->nodes,sizeof(double)); index->values=malloc(index->nodes*sizeof(double));
  if (!index->length || !index->values) return CRPX_ERROR_ALLOCATION;
  for (size_t node=0;node<index->nodes;++node) {
    index->values[node]=NAN;
    if (node!=index->root) {
      status=crpx_tree_branch_length(tree,(crpx_node_id)node,index->length+node);
      if (status!=CRPX_OK) return status;
      if (!isfinite(index->length[node]) || index->length[node]<=0) return CRPX_ERROR_INVALID_ARGUMENT;
    }
  }
  for (size_t i=0;i<count;++i) {
    if (tips[i]>=index->nodes || crpx_tree_first_child(tree,tips[i])!=CRPX_NODE_NONE ||
        !isfinite(values[i]) || isfinite(index->values[tips[i]])) return CRPX_ERROR_INVALID_ARGUMENT;
    index->values[tips[i]]=values[i];
  }
  status=crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&index->workspace);
  if (status==CRPX_OK) status=crpx_tree_workspace_refresh(index->workspace,tree);
  return status;
}
/* Gaussian-kernel product; infinite variance is flat, zero is an exact value.
 * Ratios avoid the overflow-prone product av*bv. */
static void gaussian_product(double am, double av, double bm, double bv, double *mean, double *variance)
{
  if (isinf(av)) { *mean=bm; *variance=bv; return; }
  if (isinf(bv) || av==0) { *mean=am; *variance=av; return; }
  if (bv==0) { *mean=bm; *variance=bv; return; }
  double fraction=av<=bv?(av/bv)/(1+av/bv):1/(1+bv/av);
  *mean=am*(1-fraction)+bm*fraction;
  *variance=av<=bv?av/(1+av/bv):bv/(1+bv/av);
}
void crpx_brownian_options_init(crpx_brownian_options *o) { if (o) o->reconstruct=false; }
void crpx_brownian_result_free(crpx_brownian_result *r)
{ if (r) { free(r->mean); free(r->variance); free(r); } }
size_t crpx_brownian_result_nodes(const crpx_brownian_result *r) { return r?r->nodes:0; }
size_t crpx_brownian_result_tips(const crpx_brownian_result *r) { return r?r->tips:0; }
double crpx_brownian_result_root_mean(const crpx_brownian_result *r) { return r?r->root_mean:NAN; }
double crpx_brownian_result_sigma2(const crpx_brownian_result *r) { return r?r->sigma2:NAN; }
double crpx_brownian_result_log_likelihood(const crpx_brownian_result *r) { return r?r->log_likelihood:NAN; }
const double *crpx_brownian_result_means(const crpx_brownian_result *r, size_t *count)
{ if (count) *count=r&&r->mean?r->nodes:0; return r?r->mean:NULL; }
const double *crpx_brownian_result_variances(const crpx_brownian_result *r, size_t *count)
{ if (count) *count=r&&r->variance?r->nodes:0; return r?r->variance:NULL; }
crpx_status crpx_brownian_fit(const crpx_tree *tree, const crpx_node_id *tip_nodes,
                             const double *values, size_t count, const crpx_brownian_options *options,
                             crpx_brownian_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  continuous_tree index={0}; crpx_brownian_result *r=NULL;
  double *scratch=NULL; crpx_node_id *previous=NULL,*last=NULL;
  crpx_status status=tree_prepare(tree,tip_nodes,values,count,&index);
  if (status!=CRPX_OK) goto done;
  size_t n=index.nodes; bool reconstruct=options&&options->reconstruct;
  size_t arrays=reconstruct?5:1;
  if (n>SIZE_MAX/arrays/sizeof(double) || n>SIZE_MAX/sizeof(crpx_node_id)) {
    status=CRPX_ERROR_OVERFLOW; goto done;
  }
  r=calloc(1,sizeof(*r)); scratch=calloc(n*arrays,sizeof(double));
  if (!r || !scratch) { status=CRPX_ERROR_ALLOCATION; goto done; }
  r->nodes=n; r->tips=index.tips;
  if (reconstruct) {
    previous=malloc(n*sizeof(*previous)); last=malloc(n*sizeof(*last));
    r->mean=malloc(n*sizeof(double)); r->variance=malloc(n*sizeof(double));
    if (!previous || !last || !r->mean || !r->variance) { status=CRPX_ERROR_ALLOCATION; goto done; }
  }
  double *dm=index.values,*dv=scratch;
  const crpx_node_id *postorder=crpx_tree_workspace_postorder(index.workspace);
  long double sse=0,logdet=0;
  for (size_t i=0;i<n;++i) {
    crpx_node_id node=postorder[i],child=crpx_tree_first_child(tree,node);
    if (reconstruct) last[node]=CRPX_NODE_NONE;
    if (child==CRPX_NODE_NONE) continue;
    bool first=true; crpx_node_id prev=CRPX_NODE_NONE;
    for (;child!=CRPX_NODE_NONE;child=crpx_tree_next_sibling(tree,child)) {
      if (reconstruct) { previous[child]=prev; prev=child; last[node]=child; }
      double v=dv[child]+index.length[child];
      if (!isfinite(v) || v<=0) { status=CRPX_ERROR_OVERFLOW; goto done; }
      if (first) { dm[node]=dm[child]; dv[node]=v; first=false; }
      else {
        double sum=dv[node]+v;
        long double z=((long double)dm[node]-dm[child])/sqrt(sum);
        if (!isfinite(sum) || !isfinite(z)) { status=CRPX_ERROR_OVERFLOW; goto done; }
        sse+=z*z; logdet+=log(sum);
        gaussian_product(dm[node],dv[node],dm[child],v,dm+node,dv+node);
      }
    }
  }
  logdet+=log(dv[index.root]); r->root_mean=dm[index.root];
  long double sigma2=sse/(long double)index.tips;
  if (!isfinite(sse) || !isfinite(logdet) || !isfinite(r->root_mean) || sigma2>DBL_MAX) {
    status=CRPX_ERROR_OVERFLOW; goto done;
  }
  r->sigma2=(double)sigma2;
  if (!(sigma2>0)) { status=CRPX_ERROR_INVALID_ARGUMENT; goto done; }
  if (!(r->sigma2>0)) { status=CRPX_ERROR_OVERFLOW; goto done; }
  long double log_likelihood=-0.5L*((long double)index.tips*(LOG_2PI+log(r->sigma2)+1)+logdet);
  if (!isfinite(log_likelihood) || fabsl(log_likelihood)>DBL_MAX) { status=CRPX_ERROR_OVERFLOW; goto done; }
  r->log_likelihood=(double)log_likelihood;
  if (reconstruct) {
    double *um=dv+n,*uv=um+n,*lm=uv+n,*lv=lm+n;
    uv[index.root]=INFINITY;
    const crpx_node_id *preorder=crpx_tree_workspace_preorder(index.workspace);
    for (size_t i=0;i<n;++i) {
      crpx_node_id node=preorder[i];
      gaussian_product(dm[node],dv[node],um[node],uv[node],r->mean+node,r->variance+node);
      r->variance[node]*=r->sigma2;
      if (!isfinite(r->mean[node]) || !isfinite(r->variance[node])) {
        status=CRPX_ERROR_OVERFLOW; goto done;
      }
      /* Prefix/suffix products exclude each child without subtracting nearly
       * equal precisions. Every edge participates in a bounded number of products. */
      double pm=um[node],pv=uv[node];
      for (crpx_node_id child=crpx_tree_first_child(tree,node);child!=CRPX_NODE_NONE;
           child=crpx_tree_next_sibling(tree,child)) {
        lm[child]=pm; lv[child]=pv;
        gaussian_product(pm,pv,dm[child],dv[child]+index.length[child],&pm,&pv);
      }
      double sm=0,sv=INFINITY;
      for (crpx_node_id child=last[node];child!=CRPX_NODE_NONE;child=previous[child]) {
        gaussian_product(lm[child],lv[child],sm,sv,um+child,uv+child);
        uv[child]+=index.length[child];
        gaussian_product(sm,sv,dm[child],dv[child]+index.length[child],&sm,&sv);
      }
    }
  }
  *out=r; r=NULL; status=CRPX_OK;
done:
  crpx_brownian_result_free(r); free(scratch); free(previous); free(last); tree_clear(&index); return status;
}
void crpx_pic_result_free(crpx_pic_result *r)
{ if (r) { free(r->node); free(r->contrast); free(r->variance); free(r); } }
size_t crpx_pic_result_count(const crpx_pic_result *r) { return r?r->count:0; }
const crpx_node_id *crpx_pic_result_nodes(const crpx_pic_result *r) { return r?r->node:NULL; }
const double *crpx_pic_result_contrasts(const crpx_pic_result *r) { return r?r->contrast:NULL; }
const double *crpx_pic_result_variances(const crpx_pic_result *r) { return r?r->variance:NULL; }
crpx_status crpx_pic_calculate(const crpx_tree *tree, const crpx_node_id *tip_nodes,
                              const double *values, size_t count, crpx_pic_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  continuous_tree index={0}; crpx_pic_result *r=NULL; double *variance=NULL;
  crpx_status status=tree_prepare(tree,tip_nodes,values,count,&index);
  if (status!=CRPX_OK) goto done;
  if (index.tips-1>SIZE_MAX/sizeof(crpx_node_id)) { status=CRPX_ERROR_OVERFLOW; goto done; }
  r=calloc(1,sizeof(*r)); variance=calloc(index.nodes,sizeof(double));
  if (!r || !variance) { status=CRPX_ERROR_ALLOCATION; goto done; }
  r->node=malloc((index.tips-1)*sizeof(crpx_node_id)); r->contrast=malloc((index.tips-1)*sizeof(double));
  r->variance=malloc((index.tips-1)*sizeof(double));
  if (!r->node || !r->contrast || !r->variance) { status=CRPX_ERROR_ALLOCATION; goto done; }
  const crpx_node_id *order=crpx_tree_workspace_postorder(index.workspace);
  double *mean=index.values;
  for (size_t i=0;i<index.nodes;++i) {
    crpx_node_id node=order[i],a=crpx_tree_first_child(tree,node);
    if (a==CRPX_NODE_NONE) continue;
    crpx_node_id b=crpx_tree_next_sibling(tree,a);
    if (b==CRPX_NODE_NONE || crpx_tree_next_sibling(tree,b)!=CRPX_NODE_NONE) {
      status=CRPX_ERROR_UNSUPPORTED; goto done;
    }
    double av=variance[a]+index.length[a],bv=variance[b]+index.length[b],v=av+bv;
    long double contrast=((long double)mean[a]-mean[b])/sqrt(v);
    if (!isfinite(v) || !isfinite(contrast) || fabsl(contrast)>DBL_MAX) {
      status=CRPX_ERROR_OVERFLOW; goto done;
    }
    r->node[r->count]=node; r->contrast[r->count]=(double)contrast; r->variance[r->count++]=v;
    gaussian_product(mean[a],av,mean[b],bv,mean+node,variance+node);
  }
  *out=r; r=NULL; status=CRPX_OK;
done:
  crpx_pic_result_free(r); free(variance); tree_clear(&index); return status;
}
