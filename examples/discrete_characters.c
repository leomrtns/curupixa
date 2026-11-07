/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "curupixa.h"
#include <stdio.h>

int main(void)
{
  const char *names[] = {"A", "B", "C"};
  /* Three sites, then three taxa per site; 255 is unknown. */
  const uint8_t observations[] = {0,0,1, 1,1,0, UINT8_MAX,0,1};
  crpx_tree *tree = NULL;
  crpx_character_data *data = NULL;
  crpx_substitution_model *model = NULL;
  crpx_character_analysis *analysis = NULL;
  crpx_status status = crpx_tree_from_newick("((A:0.1,B:0.2):0.3,C:0.4);", &tree, NULL);
  double log_likelihood, sitewise, unrooted_sitewise, parsimony;
  if (status == CRPX_OK) status = crpx_character_data_from_binary(names,3,3,observations,NULL,true,&data);
  if (status == CRPX_OK) status = crpx_substitution_model_poisson(2,&model);
  if (status == CRPX_OK) status = crpx_character_analysis_new(&analysis);
  /* NULL likelihood options uses the individual lengths stored on the tree. */
  if (status == CRPX_OK) status = crpx_character_analysis_likelihood(tree,data,model,NULL,analysis,&log_likelihood);
  crpx_likelihood_options options;
  crpx_likelihood_options_init(&options);
  options.branch_mode = CRPX_BRANCH_SITEWISE_EXPONENTIAL;
  options.mean_branch_length = 0.2;
  if (status == CRPX_OK) status = crpx_character_analysis_likelihood(tree,data,model,&options,analysis,&sitewise);
  /* An artificial binary root can instead represent one unrooted edge. */
  options.collapse_binary_root = true;
  options.root_edge_value = 0.2;
  if (status == CRPX_OK)
    status = crpx_character_analysis_likelihood(tree,data,model,&options,analysis,&unrooted_sitewise);
  /* NULL costs uses unordered unit changes; the tree/data/workspace stay the same. */
  if (status == CRPX_OK) status = crpx_character_analysis_parsimony(tree,data,NULL,0,analysis,&parsimony);
  if (status == CRPX_OK)
    printf("fixed log likelihood = %.9g; sitewise = %.9g; collapsed-root sitewise = %.9g; parsimony = %.9g\n",
           log_likelihood, sitewise, unrooted_sitewise, parsimony);
  else fprintf(stderr,"character analysis failed: %s\n",crpx_status_string(status));
  crpx_character_analysis_free(analysis); crpx_substitution_model_free(model);
  crpx_character_data_free(data); crpx_tree_free(tree);
  return status == CRPX_OK ? 0 : 1;
}
