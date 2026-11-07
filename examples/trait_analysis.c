/* Standalone C17/C++17 example using an installed curupixa package.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <curupixa/likelihood.h>
#include <curupixa/reconstruction.h>
#include <curupixa/continuous_traits.h>
#include <stdio.h>
#include <stdlib.h>
int main(void)
{
  crpx_tree *tree=NULL; crpx_character_data *data=NULL; crpx_substitution_model *model=NULL;
  crpx_likelihood_binding *binding=NULL; crpx_likelihood_workspace *likelihood=NULL;
  crpx_reconstruction_workspace *reconstruction=NULL; crpx_likelihood_reconstruction *asr=NULL;
  crpx_parsimony_reconstruction *parsimony=NULL; crpx_brownian_result *brownian=NULL;
  crpx_pic_result *pic=NULL; crpx_node_id tips[2];
  const char *names[]={"A","B"}; const uint8_t observations[]={0,1}; const double traits[]={1,3};
  crpx_status status=crpx_tree_from_newick("(A:0.2,B:0.3);",&tree,NULL);
  if (status!=CRPX_OK) goto done;
  status=crpx_character_data_from_binary(names,2,1,observations,NULL,false,&data);
  if (status==CRPX_OK) status=crpx_substitution_model_poisson(2,&model);
  if (status==CRPX_OK) status=crpx_likelihood_binding_new(tree,data,&binding);
  if (status==CRPX_OK) status=crpx_likelihood_workspace_new(&likelihood);
  if (status==CRPX_OK) {
    crpx_likelihood_options options; crpx_likelihood_options_init(&options); options.retain_partials=false;
    double log_likelihood;
    /* In a rate search, change options.edge_values between these evaluations.
     * The prepared binding shares data and skips repeated mapping/traversal. */
    status=crpx_likelihood_evaluate_prepared(binding,model,&options,likelihood,&log_likelihood);
    if (status==CRPX_OK) printf("log likelihood: %.12g\n",log_likelihood);
  }
  if (status==CRPX_OK) status=crpx_reconstruction_workspace_new(&reconstruction);
  if (status==CRPX_OK) status=crpx_likelihood_reconstruct(tree,data,model,NULL,reconstruction,&asr);
  if (status==CRPX_OK) {
    const double *p=crpx_likelihood_reconstruction_probabilities(asr,crpx_tree_root(tree),0);
    printf("root marginal: %.12g %.12g\n",p[0],p[1]);
    status=crpx_parsimony_reconstruct(tree,data,CRPX_PARSIMONY_UNORDERED,NULL,0,0,reconstruction,&parsimony);
  }
  if (status==CRPX_OK) {
    printf("parsimony: %.12g\n",crpx_parsimony_reconstruction_score(parsimony));
    status=crpx_character_data_map_tips(data,tree,tips,2);
  }
  if (status==CRPX_OK) status=crpx_brownian_fit(tree,tips,traits,2,NULL,&brownian); /* fit-only */
  if (status==CRPX_OK) status=crpx_pic_calculate(tree,tips,traits,2,&pic);
  if (status==CRPX_OK)
    printf("Brownian root mean: %.12g; contrasts: %zu\n",crpx_brownian_result_root_mean(brownian),
           crpx_pic_result_count(pic));
done:
  if (status!=CRPX_OK) fprintf(stderr,"curupixa: %s\n",crpx_status_string(status));
  crpx_pic_result_free(pic); crpx_brownian_result_free(brownian);
  crpx_likelihood_reconstruction_free(asr); crpx_parsimony_reconstruction_free(parsimony);
  crpx_reconstruction_workspace_free(reconstruction); crpx_likelihood_workspace_free(likelihood);
  crpx_likelihood_binding_free(binding); crpx_substitution_model_free(model);
  crpx_character_data_free(data); crpx_tree_free(tree);
  return status==CRPX_OK?EXIT_SUCCESS:EXIT_FAILURE;
}
