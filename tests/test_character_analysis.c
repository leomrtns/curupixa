/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "character_analysis.h"
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(expression) assert((expression) == CRPX_OK)

static double exhaustive(const crpx_tree *tree, const crpx_character_data *data, size_t pattern,
                         const double *costs)
{
  size_t nodes = crpx_tree_node_count(tree), states = crpx_character_data_states(data);
  size_t taxa = crpx_character_data_taxa(data), combinations = 1;
  size_t assigned[16];
  crpx_node_id tips[8];
  assert(nodes < 16 && taxa <= 8);
  OK(crpx_character_data_map_tips(data,tree,tips,8));
  for (size_t i = 0; i < nodes; ++i) combinations *= states;
  const double *values = crpx_character_data_pattern_data(data)+pattern*taxa*states;
  double minimum = INFINITY;
  for (size_t assignment = 0; assignment < combinations; ++assignment) {
    size_t code = assignment;
    for (size_t node = 0; node < nodes; ++node) { assigned[node] = code%states; code /= states; }
    bool valid = true;
    for (size_t taxon = 0; taxon < taxa; ++taxon)
      if (!values[taxon*states+assigned[tips[taxon]]]) { valid = false; break; }
    if (!valid) continue;
    double sum = 0;
    for (size_t node = 0; node < nodes; ++node) {
      crpx_node_id parent = crpx_tree_parent(tree,(crpx_node_id)node);
      if (parent != CRPX_NODE_NONE)
        sum += costs ? costs[assigned[parent]*states+assigned[node]] : assigned[parent] != assigned[node];
    }
    if (sum < minimum) minimum = sum;
  }
  return minimum;
}

static void test_exhaustive(void)
{
  const char *names[] = {"A","B","C","D"};
  size_t observations[81*4], weights[81];
  for (size_t site = 0; site < 81; ++site) {
    size_t code = site;
    for (size_t taxon = 0; taxon < 4; ++taxon) { observations[site*4+taxon] = code%3; code /= 3; }
    weights[site] = site%4;
  }
  crpx_character_data *data;
  OK(crpx_character_data_from_states(names,4,81,3,observations,weights,true,&data));
  crpx_character_analysis *analysis;
  OK(crpx_character_analysis_new(&analysis));
  const char *trees[] = {"((A,B),(C,D));","(A,B,C,D);","(((A),B),(C,D));"};
  const double costs[] = {0,2,1, 3,0,2, 1,4,0};
  for (size_t i = 0; i < 3; ++i) {
    crpx_tree *tree;
    OK(crpx_tree_from_newick(trees[i],&tree,NULL));
    for (size_t model = 0; model < 2; ++model) {
      double total, expected = 0;
      const double *matrix = model ? costs : NULL;
      OK(crpx_character_analysis_parsimony(tree,data,matrix,model ? 9 : 0,analysis,&total));
      size_t count;
      const double *scores = crpx_character_analysis_pattern_scores(analysis,&count);
      assert(count == crpx_character_data_patterns(data));
      for (size_t p = 0; p < count; ++p) {
        double best = exhaustive(tree,data,p,matrix);
        assert(scores[p] == best);
        expected += best*crpx_character_data_pattern_weights(data)[p];
      }
      assert(total == expected);
    }
    crpx_tree_free(tree); /* Reuse workspace after freeing the previous tree. */
  }
  crpx_character_analysis_free(analysis); crpx_character_data_free(data);
}

static void test_data_and_large_states(void)
{
  const char *names[] = {"A","B"};
  const size_t values[] = {69,0, 69,0, SIZE_MAX,SIZE_MAX};
  const size_t weights[] = {2,3,0};
  crpx_character_data *data;
  OK(crpx_character_data_from_states(names,2,3,70,values,weights,true,&data));
  assert(crpx_character_data_patterns(data) == 2 && crpx_character_data_total_weight(data) == 5);
  assert(crpx_character_data_site_patterns(data)[0] == crpx_character_data_site_patterns(data)[1]);
  assert(crpx_character_data_pattern_data(data)[69] == 1);
  crpx_tree *tree;
  OK(crpx_tree_from_newick("(B,A);",&tree,NULL));
  crpx_character_analysis *analysis;
  OK(crpx_character_analysis_new(&analysis));
  double score;
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score)); assert(score == 5);
  crpx_character_data_free(data);
  const uint64_t masks[] = {UINT64_C(1) << 63,1};
  OK(crpx_character_data_from_masks(names,2,1,64,masks,NULL,true,&data));
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score)); assert(score == 1);
  crpx_character_data_free(data);
  const uint8_t binary[] = {0,1, UINT8_MAX,1};
  OK(crpx_character_data_from_binary(names,2,2,binary,NULL,true,&data));
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score)); assert(score == 1);
  crpx_character_data_free(data);
  const double soft[] = {.25,.75, .5,.5};
  OK(crpx_character_data_from_tip_likelihoods(names,2,1,2,soft,NULL,true,&data));
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score)); assert(score == 0);
  crpx_character_data_free(data);
  OK(crpx_character_data_from_states(names,2,0,2,NULL,NULL,true,&data));
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score)); assert(score == 0);
  crpx_character_data_free(data); crpx_tree_free(tree); crpx_character_analysis_free(analysis);
}

static void test_failures(void)
{
  const char *names[] = {"A","B"}, *duplicates[] = {"A","A"};
  const uint8_t values[] = {0,1};
  crpx_character_data *data = (crpx_character_data *)(uintptr_t)1;
  assert(crpx_character_data_from_binary(duplicates,2,1,values,NULL,true,&data) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(!data);
  const double zeros[] = {0,0,1,0}, nan[] = {NAN,0,1,0};
  assert(crpx_character_data_from_tip_likelihoods(names,2,1,2,zeros,NULL,true,&data) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(crpx_character_data_from_tip_likelihoods(names,2,1,2,nan,NULL,true,&data) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_character_analysis *analysis;
  crpx_tree *tree;
  OK(crpx_character_analysis_new(&analysis)); OK(crpx_tree_from_newick("(A,B);",&tree,NULL));
  double score = -7;
#if SIZE_MAX > 9007199254740992ULL
  const size_t huge[] = {UINT64_C(9007199254740993)};
  OK(crpx_character_data_from_binary(names,2,1,values,huge,true,&data));
  assert(crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&score) == CRPX_ERROR_OVERFLOW && score == -7);
  size_t count = 99;
  assert(!crpx_character_analysis_pattern_scores(analysis,&count) && !count);
  crpx_character_data_free(data);
#endif
  OK(crpx_character_data_from_binary(names,2,1,values,NULL,true,&data));
  crpx_tree *bad;
  OK(crpx_tree_from_newick("(A,A);",&bad,NULL));
  crpx_node_id mapping[] = {77,88};
  assert(crpx_character_data_map_tips(data,bad,mapping,2) == CRPX_ERROR_TOPOLOGY && mapping[0] == 77 && mapping[1] == 88);
  const double invalid[] = {0,-1,1,0};
  assert(crpx_character_analysis_parsimony(tree,data,invalid,4,analysis,&score) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(bad); crpx_tree_free(tree); crpx_character_data_free(data); crpx_character_analysis_free(analysis);
}

typedef struct {
  const crpx_tree *tree;
  const crpx_character_data *data;
  const crpx_substitution_model *model;
  double expected[3];
  bool passed;
} thread_case;

static void *analyse_worker(void *argument)
{
  thread_case *test = argument;
  crpx_character_analysis *analysis = NULL;
  crpx_status status = crpx_character_analysis_new(&analysis);
  crpx_likelihood_options options;
  crpx_likelihood_options_init(&options);
  options.branch_mode = CRPX_BRANCH_SITEWISE_EXPONENTIAL;
  options.mean_branch_length = 0.17;
  options.retain_partials = false;
  for (size_t repeat = 0; repeat < 100 && status == CRPX_OK; ++repeat) {
    double scores[3];
    status = crpx_character_analysis_likelihood(test->tree,test->data,test->model,NULL,analysis,scores);
    if (status == CRPX_OK)
      status = crpx_character_analysis_likelihood(test->tree,test->data,test->model,&options,analysis,scores+1);
    if (status == CRPX_OK)
      status = crpx_character_analysis_parsimony(test->tree,test->data,NULL,0,analysis,scores+2);
    if (status == CRPX_OK) for (size_t i = 0; i < 3; ++i)
      if (scores[i] != test->expected[i]) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
  }
  test->passed = status == CRPX_OK;
  crpx_character_analysis_free(analysis);
  return NULL;
}

static void test_concurrent_analysis(void)
{
  const char *names[] = {"A","B","C","D"};
  const size_t observations[] = {0,1,2,0, 1,1,0,2, SIZE_MAX,2,0,0};
  const size_t weights[] = {1,3,2};
  crpx_tree *tree;
  crpx_character_data *data;
  crpx_substitution_model *model;
  crpx_character_analysis *serial;
  OK(crpx_tree_from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);",&tree,NULL));
  OK(crpx_character_data_from_states(names,4,3,3,observations,weights,true,&data));
  OK(crpx_substitution_model_poisson(3,&model));
  OK(crpx_character_analysis_new(&serial));
  crpx_likelihood_options options;
  crpx_likelihood_options_init(&options);
  options.branch_mode = CRPX_BRANCH_SITEWISE_EXPONENTIAL;
  options.mean_branch_length = 0.17;
  options.retain_partials = false;
  double expected[3];
  OK(crpx_character_analysis_likelihood(tree,data,model,NULL,serial,expected));
  OK(crpx_character_analysis_likelihood(tree,data,model,&options,serial,expected+1));
  OK(crpx_character_analysis_parsimony(tree,data,NULL,0,serial,expected+2));
  crpx_character_analysis_free(serial);
  pthread_t threads[4];
  thread_case cases[4];
  for (size_t worker = 0; worker < 4; ++worker) {
    cases[worker] = (thread_case){tree,data,model,{expected[0],expected[1],expected[2]},false};
    assert(!pthread_create(threads+worker,NULL,analyse_worker,cases+worker));
  }
  for (size_t worker = 0; worker < 4; ++worker) {
    assert(!pthread_join(threads[worker],NULL));
    assert(cases[worker].passed);
  }
  crpx_substitution_model_free(model); crpx_character_data_free(data); crpx_tree_free(tree);
}

int main(void)
{
  test_exhaustive(); test_data_and_large_states(); test_failures(); test_concurrent_analysis();
  puts("generic character data and parsimony: passed");
  return 0;
}
