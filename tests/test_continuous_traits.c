/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "continuous_traits.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define OK(call) assert((call)==CRPX_OK)
static void close_value(double a, long double b)
{ assert(fabsl((long double)a-b)<2e-11L*(1+fabsl(b))); }
static crpx_node_id find(const crpx_tree *t, const char *name)
{
  for (size_t i=0;i<crpx_tree_node_count(t);++i) {
    const char *s=crpx_tree_label(t,(crpx_node_id)i);
    if (s && !strcmp(s,name)) return (crpx_node_id)i;
  }
  assert(0); return CRPX_NODE_NONE;
}
/* Test-only dense pivoted inverse, independent of the pruning calculations. */
static void inverse(const long double *matrix, size_t n, long double *out, long double *logdet)
{
  long double a[32][64]={{0}}; assert(n<=32); *logdet=0;
  for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) { a[i][j]=matrix[i*n+j]; a[i][n+j]=i==j; }
  for (size_t col=0;col<n;++col) {
    size_t pivot=col;
    for (size_t i=col+1;i<n;++i) if (fabsl(a[i][col])>fabsl(a[pivot][col])) pivot=i;
    if (pivot!=col) for (size_t j=0;j<2*n;++j) {
      long double x=a[col][j]; a[col][j]=a[pivot][j]; a[pivot][j]=x;
    }
    long double value=a[col][col]; assert(value!=0); *logdet+=logl(fabsl(value));
    for (size_t j=0;j<2*n;++j) a[col][j]/=value;
    for (size_t i=0;i<n;++i) if (i!=col) {
      value=a[i][col]; for (size_t j=0;j<2*n;++j) a[i][j]-=value*a[col][j];
    }
  }
  for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) out[i*n+j]=a[i][n+j];
}
static crpx_node_id common(const crpx_tree *tree, crpx_node_id a, crpx_node_id b)
{
  for (crpx_node_id p=a;p!=CRPX_NODE_NONE;p=crpx_tree_parent(tree,p))
    for (crpx_node_id q=b;q!=CRPX_NODE_NONE;q=crpx_tree_parent(tree,q)) if (p==q) return p;
  assert(0); return CRPX_NODE_NONE;
}
static long double depth(const crpx_tree *tree, crpx_node_id node)
{
  long double result=0;
  while (node!=crpx_tree_root(tree)) {
    double length; OK(crpx_tree_branch_length(tree,node,&length)); result+=length;
    node=crpx_tree_parent(tree,node);
  }
  return result;
}
static void dense_case(const char *newick)
{
  crpx_tree *tree; OK(crpx_tree_from_newick(newick,&tree,NULL));
  size_t n=crpx_tree_node_count(tree),tips=0,internal=0; assert(n<32);
  crpx_node_id nodes[32]; size_t unknown[32]; double values[32],observed[32];
  for (size_t node=0;node<n;++node) {
    if (crpx_tree_first_child(tree,(crpx_node_id)node)==CRPX_NODE_NONE) {
      nodes[tips]=(crpx_node_id)node; values[tips]=(double)(tips*tips+tips%2)-2; observed[node]=values[tips++];
      unknown[node]=SIZE_MAX;
    } else unknown[node]=internal++;
  }
  crpx_brownian_options options; crpx_brownian_options_init(&options); assert(!options.reconstruct);
  crpx_brownian_result *score,*full;
  OK(crpx_brownian_fit(tree,nodes,values,tips,NULL,&score)); size_t count=77;
  assert(!crpx_brownian_result_means(score,&count) && !count);
  options.reconstruct=true; OK(crpx_brownian_fit(tree,nodes,values,tips,&options,&full));
  assert(crpx_brownian_result_nodes(full)==n && crpx_brownian_result_tips(full)==tips);
  close_value(crpx_brownian_result_log_likelihood(score),crpx_brownian_result_log_likelihood(full));
  /* Dense tip covariance independently estimates the fitted root mean and ML rate. */
  long double covariance[1024]={0},inv[1024],logdet;
  for (size_t i=0;i<tips;++i) for (size_t j=0;j<tips;++j)
    covariance[i*tips+j]=depth(tree,common(tree,nodes[i],nodes[j]));
  inverse(covariance,tips,inv,&logdet);
  long double numerator=0,denominator=0;
  for (size_t i=0;i<tips;++i) for (size_t j=0;j<tips;++j) {
    numerator+=inv[i*tips+j]*values[j]; denominator+=inv[i*tips+j];
  }
  long double mu=numerator/denominator,sse=0;
  for (size_t i=0;i<tips;++i) for (size_t j=0;j<tips;++j)
    sse+=(values[i]-mu)*inv[i*tips+j]*(values[j]-mu);
  long double sigma2=sse/tips,ll=-.5L*(tips*(logl(2*acosl(-1.0L))+logl(sigma2)+1)+logdet);
  close_value(crpx_brownian_result_root_mean(full),mu);
  close_value(crpx_brownian_result_sigma2(full),sigma2);
  close_value(crpx_brownian_result_log_likelihood(full),ll);
  /* Conditional internal-node precision from independent edge energies, using
   * a flat prior at the root; invert to obtain all marginal means/variances. */
  long double precision[1024]={0},rhs[32]={0};
  for (size_t node=0;node<n;++node) if (node!=crpx_tree_root(tree)) {
    crpx_node_id parent=crpx_tree_parent(tree,(crpx_node_id)node); double length;
    OK(crpx_tree_branch_length(tree,(crpx_node_id)node,&length)); long double rate=1.0L/length;
    size_t a=unknown[node],b=unknown[parent];
    if (a!=SIZE_MAX) precision[a*internal+a]+=rate;
    if (b!=SIZE_MAX) precision[b*internal+b]+=rate;
    if (a!=SIZE_MAX && b!=SIZE_MAX) precision[a*internal+b]-=rate,precision[b*internal+a]-=rate;
    else if (a!=SIZE_MAX) rhs[a]+=rate*observed[parent];
    else if (b!=SIZE_MAX) rhs[b]+=rate*observed[node];
  }
  inverse(precision,internal,inv,&logdet);
  const double *means=crpx_brownian_result_means(full,&count);
  assert(means && count==n); const double *variance=crpx_brownian_result_variances(full,&count);
  for (size_t node=0;node<n;++node) {
    size_t a=unknown[node];
    if (a==SIZE_MAX) { close_value(means[node],observed[node]); assert(variance[node]==0); }
    else {
      long double expected=0; for (size_t b=0;b<internal;++b) expected+=inv[a*internal+b]*rhs[b];
      close_value(means[node],expected); close_value(variance[node],inv[a*internal+a]*sigma2);
    }
  }
  uint64_t revision=crpx_tree_revision(tree);
  crpx_pic_result *pic=NULL;
  crpx_status status=crpx_pic_calculate(tree,nodes,values,tips,&pic);
  assert(status==CRPX_OK || status==CRPX_ERROR_UNSUPPORTED);
  if (status==CRPX_OK) {
    long double sum=0; const double *contrasts=crpx_pic_result_contrasts(pic);
    for (size_t i=0;i<crpx_pic_result_count(pic);++i) sum+=(long double)contrasts[i]*contrasts[i];
    assert(crpx_pic_result_count(pic)==tips-1); close_value((double)sum,sse);
  }
  assert(crpx_tree_revision(tree)==revision);
  crpx_tree_free(tree); assert(isfinite(means[0])); /* result owns all arrays */
  crpx_brownian_result_free(score); crpx_brownian_result_free(full); crpx_pic_result_free(pic);
}
static void test_pic_and_failures(void)
{
  crpx_tree *tree; OK(crpx_tree_from_newick("((A:1,B:3):2,C:5):-7;",&tree,NULL));
  crpx_node_id nodes[]={find(tree,"C"),find(tree,"A"),find(tree,"B")}; double values[]={8,1,5};
  crpx_pic_result *pic; OK(crpx_pic_calculate(tree,nodes,values,3,&pic));
  const double *contrasts=crpx_pic_result_contrasts(pic),*variance=crpx_pic_result_variances(pic);
  close_value(contrasts[0],-2); close_value(variance[0],4);
  close_value(contrasts[1],-6/sqrtl(7.75L)); close_value(variance[1],7.75L);
  crpx_pic_result_free(pic);
  crpx_brownian_result *r=(void *)1; crpx_node_id duplicate[]={nodes[0],nodes[0],nodes[2]};
  assert(crpx_brownian_fit(tree,duplicate,values,3,NULL,&r)==CRPX_ERROR_INVALID_ARGUMENT && !r);
  values[0]=NAN;
  assert(crpx_brownian_fit(tree,nodes,values,3,NULL,&r)==CRPX_ERROR_INVALID_ARGUMENT && !r);
  values[0]=values[1]=values[2]=4;
  assert(crpx_brownian_fit(tree,nodes,values,3,NULL,&r)==CRPX_ERROR_INVALID_ARGUMENT && !r);
  OK(crpx_pic_calculate(tree,nodes,values,3,&pic));
  for (size_t i=0;i<crpx_pic_result_count(pic);++i) assert(crpx_pic_result_contrasts(pic)[i]==0);
  crpx_pic_result_free(pic);
  OK(crpx_tree_set_branch_length(tree,nodes[0],0));
  assert(crpx_pic_calculate(tree,nodes,values,3,&pic)==CRPX_ERROR_INVALID_ARGUMENT && !pic);
  OK(crpx_tree_clear_branch_length(tree,nodes[0]));
  assert(crpx_pic_calculate(tree,nodes,values,3,&pic)==CRPX_ERROR_NOT_FOUND && !pic);
  crpx_tree_free(tree);
}
int main(void)
{
  dense_case("((A:1,B:3):2,C:5);");
  dense_case("(A:1,B:2,C:3,D:4);");
  dense_case("(((A:1,B:2):3):4,(C:2,D:5):1);");
  dense_case("(((A:1,B:2):3,C:2):1,D:5):-999;");
  test_pic_and_failures();
  puts("linear-memory Brownian/PIC vs independent dense Gaussian references: passed"); return 0;
}
