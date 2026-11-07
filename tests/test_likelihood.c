/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "likelihood.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(call) assert((call) == CRPX_OK)
static void close_value(double a, double b) { assert(fabs(a-b) < 1e-10*(1+fabs(b))); }

/* Independent full enumeration over every node state (not a pruning recursion). */
static double enumerate(const crpx_tree *tree, size_t states, const double *pi, const double *matrices,
                         const double *tips, const size_t *tip_index)
{
  size_t nodes = crpx_tree_node_count(tree), assignments = 1;
  for (size_t i = 0; i < nodes; ++i) assignments *= states;
  size_t *assignment = malloc(nodes*sizeof(size_t)); assert(assignment);
  long double total = 0;
  for (size_t code = 0; code < assignments; ++code) {
    size_t value = code;
    for (size_t i = 0; i < nodes; ++i) { assignment[i] = value%states; value /= states; }
    long double probability = pi[assignment[crpx_tree_root(tree)]];
    for (size_t node = 0; node < nodes; ++node) {
      crpx_node_id parent = crpx_tree_parent(tree,(crpx_node_id)node);
      if (parent != CRPX_NODE_NONE)
        probability *= matrices[(node*states+assignment[parent])*states+assignment[node]];
      if (tip_index[node] != SIZE_MAX) probability *= tips[tip_index[node]*states+assignment[node]];
    }
    total += probability;
  }
  free(assignment); return (double)total;
}
static void exhaustive_case(const char *newick, size_t states, bool hky, crpx_branch_mode mode)
{
  crpx_tree *tree; crpx_character_data *data; crpx_substitution_model *model;
  crpx_likelihood_workspace *w;
  OK(crpx_tree_from_newick(newick,&tree,NULL));
  const char *names[] = {"A","B","C"};
  size_t nodes = crpx_tree_node_count(tree), observations[27], weight[] = {1,2,3,4,5,6,7,8,9};
  for (size_t site = 0; site < 9; ++site) for (size_t taxon = 0; taxon < 3; ++taxon)
    observations[site*3+taxon] = (site+site/3+taxon*2)%states;
  observations[3] = SIZE_MAX;
  OK(crpx_character_data_from_states(names,3,9,states,observations,weight,false,&data));
  double frequencies[] = {.13,.29,.31,.27};
  if (hky) OK(crpx_substitution_model_hky(frequencies,2.7,&model));
  else OK(crpx_substitution_model_poisson(states,&model));
  OK(crpx_likelihood_workspace_new(&w));
  double *matrices = calloc(nodes*states*states,sizeof(double)), *edges = malloc(nodes*sizeof(double));
  size_t *tip_index = malloc(nodes*sizeof(size_t)); crpx_node_id mapped[3];
  assert(matrices && edges && tip_index);
  OK(crpx_character_data_map_tips(data,tree,mapped,3));
  for (size_t i = 0; i < nodes; ++i) { tip_index[i] = SIZE_MAX; edges[i] = .1+.03*i; }
  for (size_t i = 0; i < 3; ++i) tip_index[mapped[i]] = i;
  for (size_t i = 0; i < nodes; ++i)
    OK(crpx_substitution_model_transition(model,edges[i],mode,matrices+i*states*states,states*states));
  crpx_likelihood_options o; crpx_likelihood_options_init(&o);
  o.branch_mode = mode; o.edge_values = edges; o.edge_count = nodes;
  edges[crpx_tree_root(tree)] = NAN; /* stem ignored */
  double total; OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&total));
  const double *patterns = crpx_character_data_pattern_data(data);
  const double *logs = crpx_likelihood_workspace_pattern_log_likelihoods(w,NULL);
  double reference = 0;
  for (size_t p = 0; p < 9; ++p) {
    double probability = enumerate(tree,states,crpx_substitution_model_frequencies(model),matrices,
                                   patterns+p*3*states,tip_index);
    close_value(logs[p],log(probability)); reference += weight[p]*log(probability);
  }
  close_value(total,reference);
  const double *partial; double scale;
  OK(crpx_likelihood_workspace_partial(w,crpx_tree_root(tree),0,&partial,&scale));
  double root_sum = 0;
  for (size_t s = 0; s < states; ++s) root_sum += crpx_substitution_model_frequencies(model)[s]*partial[s];
  close_value(log(root_sum)+scale,logs[0]);
  o.retain_partials = false;
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&total)); close_value(total,reference);
  assert(crpx_likelihood_workspace_partial(w,0,0,&partial,&scale) == CRPX_ERROR_NOT_FOUND);
  edges[mapped[0]] = -1; total = 123;
  assert(crpx_likelihood_evaluate(tree,data,model,&o,w,&total) == CRPX_ERROR_INVALID_ARGUMENT && total == 123);
  size_t count = 999;
  assert(!crpx_likelihood_workspace_pattern_log_likelihoods(w,&count) && count == 0);
  free(matrices); free(edges); free(tip_index); crpx_likelihood_workspace_free(w);
  crpx_character_data_free(data); crpx_substitution_model_free(model); crpx_tree_free(tree);
}
static void test_reference_semantics(void)
{
  crpx_tree *tree; crpx_character_data *data; crpx_substitution_model *model; crpx_likelihood_workspace *w;
  const char *names[] = {"A","B"}; uint8_t observations[] = {0,0,0,0};
  OK(crpx_tree_from_newick("(A,B);",&tree,NULL));
  OK(crpx_character_data_from_binary(names,2,2,observations,NULL,true,&data));
  OK(crpx_substitution_model_poisson(2,&model)); OK(crpx_likelihood_workspace_new(&w));
  crpx_likelihood_options o; crpx_likelihood_options_init(&o);
  o.branch_mode = CRPX_BRANCH_SITEWISE_EXPONENTIAL; o.mean_branch_length = .5;
  o.collapse_binary_root = true; o.root_edge_value = .5;
  double score; OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score));
  close_value(exp(score),9.0/64.0); assert(fabs(exp(score)-7.0/48.0) > .001);
  o.collapse_binary_root = false;
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score)); close_value(exp(score),25.0/256.0);
  o.branch_mode = CRPX_BRANCH_FIXED; o.collapse_binary_root = true;
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score));
  close_value(score,2*log((1+exp(-1.0))/4));
  crpx_likelihood_workspace_free(w); crpx_substitution_model_free(model); crpx_character_data_free(data);
  uint8_t conflict[] = {0,1}; size_t zero = 0;
  OK(crpx_character_data_from_binary(names,2,1,conflict,&zero,true,&data));
  OK(crpx_substitution_model_poisson(2,&model)); OK(crpx_likelihood_workspace_new(&w));
  o.root_edge_value = 0;
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score)); assert(score == 0);
  assert(crpx_likelihood_workspace_pattern_log_likelihoods(w,NULL)[0] == -INFINITY);
  crpx_character_data_free(data);
  OK(crpx_character_data_from_binary(names,2,1,conflict,NULL,true,&data));
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score)); assert(score == -INFINITY);
  crpx_character_data_free(data); crpx_substitution_model_free(model);
  const double pi[] = {.25,.25,.25,.25}; size_t dna[] = {0,1};
  OK(crpx_character_data_from_states(names,2,1,4,dna,NULL,true,&data));
  OK(crpx_substitution_model_hky(pi,1,&model));
  o.root_edge_value = DBL_TRUE_MIN;
  OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score));
  assert(isfinite(score)); close_value(score,log(.25/3)+log(DBL_TRUE_MIN));
  crpx_tree_free(tree); /* workspace owns partials, no borrowed tree */
  const double *partial; double scale;
  OK(crpx_likelihood_workspace_log_partial(w,0,0,&partial,&scale));
  crpx_likelihood_workspace_free(w); crpx_character_data_free(data); crpx_substitution_model_free(model);
}
static void test_deep_probability(void)
{
  const size_t taxa = 2000;
  crpx_tree *tree; crpx_node_id root, node;
  char **names = calloc(taxa,sizeof(*names)); uint8_t *states = malloc(taxa); assert(names && states);
  OK(crpx_tree_new(&tree)); OK(crpx_tree_add_node(tree,CRPX_NODE_NONE,&root));
  for (size_t i = 0; i < taxa; ++i) {
    names[i] = malloc(32); assert(names[i]); snprintf(names[i],32,"taxon%zu",i); states[i] = i%2;
    OK(crpx_tree_add_node(tree,root,&node)); OK(crpx_tree_set_label(tree,node,names[i]));
    OK(crpx_tree_set_branch_length(tree,node,1e-200));
  }
  crpx_character_data *data; crpx_substitution_model *model; crpx_likelihood_workspace *w;
  OK(crpx_character_data_from_binary((const char *const *)names,taxa,1,states,NULL,true,&data));
  OK(crpx_substitution_model_poisson(2,&model)); OK(crpx_likelihood_workspace_new(&w));
  double score; OK(crpx_likelihood_evaluate(tree,data,model,NULL,w,&score));
  close_value(score,(double)(taxa/2)*log(1e-200));
  crpx_likelihood_workspace_free(w); crpx_character_data_free(data); crpx_substitution_model_free(model);
  crpx_tree_free(tree); for (size_t i = 0; i < taxa; ++i) free(names[i]); free(names); free(states);
}
int main(void)
{
  for (unsigned mode = 0; mode < 2; ++mode) {
    exhaustive_case("((A:0.1,B:0.2):0.3,C:0.4);",2,false,(crpx_branch_mode)mode);
    exhaustive_case("(A,B,C);",3,false,(crpx_branch_mode)mode);
    exhaustive_case("((A),B,C);",4,true,(crpx_branch_mode)mode);
  }
  test_reference_semantics(); test_deep_probability();
  puts("likelihood enumeration, reference semantics, partials and log-space scaling: passed");
  return 0;
}
