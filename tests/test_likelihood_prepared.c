/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "likelihood.h"
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define OK(x) assert((x)==CRPX_OK)
#ifdef CRPX_TEST_ALLOCATION_WRAPPERS
static size_t allocations,mappings,validations;
static bool counting;
void *__real_malloc(size_t); void *__real_calloc(size_t,size_t); void *__real_realloc(void *,size_t);
void *__wrap_malloc(size_t n) { if(counting) ++allocations; return __real_malloc(n); }
void *__wrap_calloc(size_t n,size_t s) { if(counting) ++allocations; return __real_calloc(n,s); }
void *__wrap_realloc(void *p,size_t n) { if(counting) ++allocations; return __real_realloc(p,n); }
crpx_status __real_crpx_character_data_map_tips(const crpx_character_data *,const crpx_tree *,crpx_node_id *,size_t);
crpx_status __wrap_crpx_character_data_map_tips(const crpx_character_data *d,const crpx_tree *t,crpx_node_id *ids,size_t n)
{ if(counting) ++mappings; return __real_crpx_character_data_map_tips(d,t,ids,n); }
crpx_status __real_crpx_tree_validate(const crpx_tree *);
crpx_status __wrap_crpx_tree_validate(const crpx_tree *t)
{ if(counting) ++validations; return __real_crpx_tree_validate(t); }
#endif
static void close_value(double a,double b) { assert(fabs(a-b)<1e-12*(1+fabs(b))); }
typedef struct { const crpx_likelihood_binding *binding; const crpx_substitution_model *model; } worker_input;
static void *worker(void *argument)
{
  const worker_input *input=argument; crpx_likelihood_workspace *w; OK(crpx_likelihood_workspace_new(&w));
  for(size_t i=0;i<30;++i) {
    double score; OK(crpx_likelihood_evaluate_prepared(input->binding,input->model,NULL,w,&score));
    close_value(score,log(-expm1(-1.0)/4));
  }
  crpx_likelihood_workspace_free(w); return NULL;
}
int main(void)
{
  const char *names[]={"A","B"}; uint8_t observations[]={0,1};
  crpx_tree *tree; crpx_character_data *data; crpx_substitution_model *model,*other;
  crpx_likelihood_binding *binding,*second; crpx_likelihood_workspace *w;
  OK(crpx_tree_from_newick("(A:0.2,B:0.3);",&tree,NULL));
  OK(crpx_character_data_from_binary(names,2,1,observations,NULL,false,&data));
  OK(crpx_substitution_model_poisson(2,&model)); double frequencies[]={.3,.7};
  OK(crpx_substitution_model_f81(2,frequencies,&other));
  OK(crpx_likelihood_workspace_new(&w)); size_t empty=crpx_likelihood_workspace_storage_size(w);
  OK(crpx_likelihood_binding_new(tree,data,&binding));
  crpx_likelihood_options o; crpx_likelihood_options_init(&o); o.retain_partials=false;
  double score,reference; OK(crpx_likelihood_evaluate_prepared(binding,model,&o,w,&score));
  close_value(score,log(-expm1(-1.0)/4));
  double edges[]={NAN,.2,.3}; o.edge_values=edges; o.edge_count=3;
#ifdef CRPX_TEST_ALLOCATION_WRAPPERS
  counting=true;
#endif
  for(size_t i=1;i<=100;++i) {
    edges[1]=.002*i; edges[2]=.003*i;
    OK(crpx_likelihood_evaluate_prepared(binding,model,&o,w,&score));
    close_value(score,log(-expm1(-2*(edges[1]+edges[2]))/4));
  }
#ifdef CRPX_TEST_ALLOCATION_WRAPPERS
  counting=false; assert(!allocations && !mappings && !validations);
  counting=true;
  for(size_t i=0;i<100;++i) OK(crpx_likelihood_evaluate(tree,data,model,&o,w,&score));
  counting=false;
  assert(allocations && mappings && validations);
  printf("100 rate evaluations: prepared 0 allocations/maps/validations; ordinary %zu/%zu/%zu\n",
         allocations,mappings,validations);
#endif
  OK(crpx_likelihood_evaluate(tree,data,other,&o,w,&reference));
  OK(crpx_likelihood_evaluate_prepared(binding,other,&o,w,&score)); close_value(score,reference);
  OK(crpx_tree_set_branch_length(tree,1,9)); OK(crpx_tree_set_label(tree,1,"changed"));
  crpx_likelihood_binding *invalid=(void *)1;
  assert(crpx_likelihood_binding_new(tree,data,&invalid)==CRPX_ERROR_TOPOLOGY && !invalid);
  OK(crpx_tree_set_label(tree,1,"A")); OK(crpx_likelihood_binding_new(tree,data,&second));
  OK(crpx_likelihood_evaluate_prepared(second,model,NULL,w,&reference));
  crpx_tree_free(tree); crpx_character_data_free(data); /* bindings own safe independent lifetimes */
  OK(crpx_likelihood_evaluate_prepared(second,model,NULL,w,&score)); close_value(score,reference);
  OK(crpx_likelihood_evaluate_prepared(binding,model,NULL,w,&score)); close_value(score,log(-expm1(-1.0)/4));
  pthread_t threads[3]; worker_input input={binding,model};
  for(size_t i=0;i<3;++i) assert(!pthread_create(threads+i,NULL,worker,&input));
  for(size_t i=0;i<3;++i) assert(!pthread_join(threads[i],NULL));
  o.edge_count=2; score=123;
  assert(crpx_likelihood_evaluate_prepared(binding,model,&o,w,&score)==CRPX_ERROR_INVALID_ARGUMENT && score==123);
  assert(!crpx_likelihood_workspace_pattern_log_likelihoods(w,NULL));
  assert(crpx_likelihood_workspace_storage_size(w)>empty); crpx_likelihood_workspace_clear(w);
  assert(crpx_likelihood_workspace_storage_size(w)==empty);
  OK(crpx_likelihood_evaluate_prepared(binding,model,NULL,w,&score));
  crpx_likelihood_binding_free(binding); crpx_likelihood_binding_free(second);
  crpx_likelihood_workspace_free(w); crpx_substitution_model_free(model); crpx_substitution_model_free(other);
  puts("prepared likelihood: analytic scores, allocation-free repeats, snapshot ownership and shared workers passed");
  return 0;
}
