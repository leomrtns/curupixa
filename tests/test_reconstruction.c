/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "reconstruction.h"
#include "character_analysis.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(x) assert((x)==CRPX_OK)
static bool close_to(double a,double b) { return fabs(a-b)<2e-10*fmax(1,fmax(fabs(a),fabs(b))); }
static size_t power(size_t k,size_t n) { size_t r=1; while (n--) r*=k; return r; }
static void states_of(size_t code,size_t *states,size_t n,size_t k)
{ for (size_t i=0;i<n;++i) { states[i]=code%k; code/=k; } }
static crpx_character_data *data_for(const crpx_tree *t,size_t k)
{
  size_t taxa=crpx_tree_leaf_count(t),j=0;
  const char *names[8]; double values[2*8*4]; size_t weights[2]={2,3};
  assert(taxa<=8 && k<=4);
  for (size_t i=0;i<crpx_tree_node_count(t);++i)
    if (crpx_tree_first_child(t,(crpx_node_id)i)==CRPX_NODE_NONE) names[j++]=crpx_tree_label(t,(crpx_node_id)i);
  for (size_t p=0;p<2;++p) for (size_t taxon=0;taxon<taxa;++taxon) for (size_t s=0;s<k;++s)
    values[(p*taxa+taxon)*k+s]=taxon+1==taxa?(p?(double)(s+1)/4:1):((taxon+p)%k==s?1:0);
  crpx_character_data *data;
  OK(crpx_character_data_from_tip_likelihoods(names,taxa,2,k,values,weights,false,&data)); return data;
}
static void likelihood_case(crpx_reconstruction_workspace *w,const char *newick,
                             const crpx_substitution_model *model,bool exponential,bool collapse)
{
  crpx_tree *t; OK(crpx_tree_from_newick(newick,&t,NULL));
  size_t n=crpx_tree_node_count(t),k=crpx_substitution_model_states(model),taxa=crpx_tree_leaf_count(t);
  assert(n<=8 && k<=4);
  crpx_character_data *data=data_for(t,k);
  crpx_likelihood_options o; crpx_likelihood_options_init(&o);
  o.branch_mode=exponential?CRPX_BRANCH_SITEWISE_EXPONENTIAL:CRPX_BRANCH_FIXED;
  o.mean_branch_length=.37; o.root_edge_value=.63; o.collapse_binary_root=collapse; o.retain_partials=false;
  crpx_likelihood_reconstruction *r;
  OK(crpx_likelihood_reconstruct(t,data,model,&o,w,&r));
  assert(crpx_likelihood_reconstruction_nodes(r)==n && crpx_likelihood_reconstruction_states(r)==k &&
    crpx_likelihood_reconstruction_patterns(r)==2);
  double kernels[8*16]; size_t assignment[8]; crpx_node_id tips[8],root=crpx_tree_root(t);
  OK(crpx_character_data_map_tips(data,t,tips,taxa));
  crpx_node_id first=CRPX_NODE_NONE,second=CRPX_NODE_NONE;
  if (collapse) { first=crpx_tree_first_child(t,root); second=crpx_tree_next_sibling(t,first); }
  for (size_t node=0;node<n;++node) if (node!=root) {
    double edge;
    if (node==first) edge=o.root_edge_value;
    else if (node==second) edge=0;
    else if (exponential) edge=o.mean_branch_length;
    else OK(crpx_tree_branch_length(t,(crpx_node_id)node,&edge));
    OK(crpx_substitution_model_transition(model,edge,o.branch_mode,kernels+node*k*k,k*k));
  }
  long double expected_score=0;
  const double *obs=crpx_character_data_pattern_data(data),*pi=crpx_substitution_model_frequencies(model);
  const size_t *weights=crpx_character_data_pattern_weights(data);
  for (size_t p=0;p<2;++p) {
    long double sum=0,marginal[8*4]={0};
    for (size_t code=0;code<power(k,n);++code) {
      states_of(code,assignment,n,k); long double probability=pi[assignment[root]];
      for (size_t node=0;node<n;++node) if (node!=root) {
        crpx_node_id parent=crpx_tree_parent(t,(crpx_node_id)node);
        probability*=kernels[node*k*k+assignment[parent]*k+assignment[node]];
      }
      for (size_t taxon=0;taxon<taxa;++taxon) probability*=obs[(p*taxa+taxon)*k+assignment[tips[taxon]]];
      sum+=probability;
      for (size_t node=0;node<n;++node) marginal[node*k+assignment[node]]+=probability;
    }
    assert(sum>0); expected_score+=weights[p]*logl(sum);
    for (size_t node=0;node<n;++node) {
      const double *prob=crpx_likelihood_reconstruction_probabilities(r,(crpx_node_id)node,p);
      double total=0;
      for (size_t s=0;s<k;++s) { assert(close_to(prob[s],(double)(marginal[node*k+s]/sum))); total+=prob[s]; }
      assert(close_to(total,1));
    }
  }
  assert(close_to(crpx_likelihood_reconstruction_score(r),(double)expected_score));
  crpx_likelihood_workspace *score_workspace; double score;
  OK(crpx_likelihood_workspace_new(&score_workspace)); OK(crpx_likelihood_evaluate(t,data,model,&o,score_workspace,&score));
  assert(close_to(score,crpx_likelihood_reconstruction_score(r)));
  crpx_likelihood_workspace_free(score_workspace); crpx_character_data_free(data); crpx_tree_free(t);
  crpx_reconstruction_workspace_clear(w);
  assert(isfinite(crpx_likelihood_reconstruction_probabilities(r,0,0)[0])); /* independent lifetime */
  crpx_likelihood_reconstruction_free(r);
}
static double cost(size_t a,size_t b,crpx_parsimony_cost_mode mode,const double *matrix,size_t k)
{ return mode==CRPX_PARSIMONY_CUSTOM?matrix[a*k+b]:mode==CRPX_PARSIMONY_ORDERED?(double)(a>b?a-b:b-a):(a!=b); }
static void parsimony_case(crpx_reconstruction_workspace *w,const char *newick,size_t k,
                            crpx_parsimony_cost_mode mode,const double *matrix)
{
  crpx_tree *t; OK(crpx_tree_from_newick(newick,&t,NULL));
  size_t n=crpx_tree_node_count(t),taxa=crpx_tree_leaf_count(t),assignment[8]; assert(n<=8);
  crpx_character_data *data=data_for(t,k); crpx_node_id tips[8],root=crpx_tree_root(t);
  OK(crpx_character_data_map_tips(data,t,tips,taxa));
  crpx_parsimony_reconstruction *r; OK(crpx_parsimony_reconstruct(t,data,mode,matrix,matrix?k*k:0,0,w,&r));
  const double *obs=crpx_character_data_pattern_data(data); const size_t *weights=crpx_character_data_pattern_weights(data);
  double total=0;
  for (size_t p=0;p<2;++p) {
    double best=INFINITY; unsigned char admissible[8*4]={0};
    for (size_t code=0;code<power(k,n);++code) {
      states_of(code,assignment,n,k); bool allowed=true;
      for (size_t taxon=0;taxon<taxa;++taxon)
        if (!obs[(p*taxa+taxon)*k+assignment[tips[taxon]]]) allowed=false;
      if (!allowed) continue;
      double candidate=0;
      for (size_t node=0;node<n;++node) if (node!=root)
        candidate+=cost(assignment[crpx_tree_parent(t,(crpx_node_id)node)],assignment[node],mode,matrix,k);
      if (candidate<best) { best=candidate; memset(admissible,0,sizeof(admissible)); }
      if (candidate==best) for (size_t node=0;node<n;++node) admissible[node*k+assignment[node]]=1;
    }
    assert(crpx_parsimony_reconstruction_pattern_scores(r)[p]==best); total+=weights[p]*best;
    double joint=0;
    for (size_t node=0;node<n;++node) {
      const unsigned char *actual=crpx_parsimony_reconstruction_optimal_states(r,(crpx_node_id)node,p);
      assert(!memcmp(actual,admissible+node*k,k)); assignment[node]=crpx_parsimony_reconstruction_assignment(r,(crpx_node_id)node,p);
      assert(assignment[node]<k && actual[assignment[node]]);
    }
    for (size_t node=0;node<n;++node) if (node!=root)
      joint+=cost(assignment[crpx_tree_parent(t,(crpx_node_id)node)],assignment[node],mode,matrix,k);
    for (size_t taxon=0;taxon<taxa;++taxon) assert(obs[(p*taxa+taxon)*k+assignment[tips[taxon]]]>0);
    assert(joint==best);
  }
  assert(crpx_parsimony_reconstruction_score(r)==total);
  crpx_character_analysis *analysis; OK(crpx_character_analysis_new(&analysis));
  double explicit_cost[16],score;
  for (size_t a=0;a<k;++a) for (size_t b=0;b<k;++b) explicit_cost[a*k+b]=cost(a,b,mode,matrix,k);
  OK(crpx_character_analysis_parsimony(t,data,mode==CRPX_PARSIMONY_UNORDERED?NULL:explicit_cost,
      mode==CRPX_PARSIMONY_UNORDERED?0:k*k,analysis,&score)); assert(score==total);
  crpx_character_analysis_free(analysis); crpx_character_data_free(data); crpx_tree_free(t);
  crpx_reconstruction_workspace_clear(w);
  assert(crpx_parsimony_reconstruction_assignment(r,0,0)<k); crpx_parsimony_reconstruction_free(r);
}
static void test_errors_and_tiny(crpx_reconstruction_workspace *w)
{
  crpx_tree *t; OK(crpx_tree_from_newick("(A:0,B:5e-324);",&t,NULL));
  const char *names[]={"A","B"}; size_t observations[]={0,1}; crpx_character_data *data;
  OK(crpx_character_data_from_states(names,2,1,2,observations,NULL,false,&data));
  crpx_substitution_model *m; OK(crpx_substitution_model_poisson(2,&m)); crpx_likelihood_reconstruction *r;
  OK(crpx_likelihood_reconstruct(t,data,m,NULL,w,&r)); assert(isfinite(crpx_likelihood_reconstruction_score(r)));
  assert(crpx_likelihood_reconstruction_probabilities(r,crpx_tree_root(t),0)[0]==1);
  crpx_likelihood_reconstruction_free(r);
  OK(crpx_tree_set_branch_length(t,2,0)); r=(void *)1;
  assert(crpx_likelihood_reconstruct(t,data,m,NULL,w,&r)==CRPX_ERROR_INVALID_ARGUMENT && !r);
  crpx_parsimony_reconstruction *pr=(void *)1;
  double bad[]={0,-1,1,0};
  assert(crpx_parsimony_reconstruct(t,data,CRPX_PARSIMONY_CUSTOM,bad,4,0,w,&pr)==CRPX_ERROR_INVALID_ARGUMENT && !pr);
  crpx_substitution_model_free(m); crpx_character_data_free(data); crpx_tree_free(t);
}
static void test_arbitrary_states(crpx_reconstruction_workspace *w,size_t k)
{
  crpx_tree *t; OK(crpx_tree_from_newick("(A:.2,B:.8);",&t,NULL));
  const char *names[]={"A","B"}; size_t observations[]={0,k-1};
  crpx_character_data *data; crpx_substitution_model *model;
  OK(crpx_character_data_from_states(names,2,1,k,observations,NULL,false,&data));
  OK(crpx_substitution_model_poisson(k,&model));
  crpx_likelihood_reconstruction *r; OK(crpx_likelihood_reconstruct(t,data,model,NULL,w,&r));
  const double *prob=crpx_likelihood_reconstruction_probabilities(r,crpx_tree_root(t),0);
  double sum=0; for (size_t s=0;s<k;++s) { assert(isfinite(prob[s]) && prob[s]>=0); sum+=prob[s]; }
  assert(close_to(sum,1)); crpx_likelihood_reconstruction_free(r);
  crpx_parsimony_reconstruction *pr;
  OK(crpx_parsimony_reconstruct(t,data,CRPX_PARSIMONY_ORDERED,NULL,0,0,w,&pr));
  assert(crpx_parsimony_reconstruction_score(pr)==(double)(k-1));
  const unsigned char *optimal=crpx_parsimony_reconstruction_optimal_states(pr,crpx_tree_root(t),0);
  for (size_t s=0;s<k;++s) assert(optimal[s]);
  crpx_parsimony_reconstruction_free(pr);
  crpx_substitution_model_free(model); crpx_character_data_free(data); crpx_tree_free(t);
}
int main(void)
{
  crpx_reconstruction_workspace *w; OK(crpx_reconstruction_workspace_new(&w));
  crpx_substitution_model *er,*f81,*hky;
  double frequencies[]={.1,.2,.3,.4};
  OK(crpx_substitution_model_poisson(3,&er)); OK(crpx_substitution_model_f81(3,frequencies,&f81));
  OK(crpx_substitution_model_hky(frequencies,2.7,&hky));
  const char *trees[]={"((A:0,B:.3):.2,C:.7);","(A:0,B:.3,C:.7);","((A:.2):.1,B:.8);","A;"};
  const crpx_substitution_model *models[]={er,f81,hky};
  for (size_t i=0;i<4;++i) for (size_t m=0;m<3;++m) {
    likelihood_case(w,trees[i],models[m],false,false); likelihood_case(w,trees[i],models[m],true,false);
    if (i==0 || i==2) likelihood_case(w,trees[i],models[m],false,true);
  }
  const double directed[]={0,2,3,1,0,4,6,1,0};
  for (size_t i=0;i<4;++i) {
    parsimony_case(w,trees[i],3,CRPX_PARSIMONY_UNORDERED,NULL);
    parsimony_case(w,trees[i],3,CRPX_PARSIMONY_ORDERED,NULL);
    parsimony_case(w,trees[i],3,CRPX_PARSIMONY_CUSTOM,directed);
  }
  test_errors_and_tiny(w);
  test_arbitrary_states(w,1); test_arbitrary_states(w,70);
  crpx_substitution_model_free(er); crpx_substitution_model_free(f81); crpx_substitution_model_free(hky);
  crpx_reconstruction_workspace_free(w);
  puts("Full-data likelihood marginals, global parsimony sets and coherent assignments: passed"); return 0;
}
