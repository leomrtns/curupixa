/* Compile against an installed package:
 * cc -std=c17 gene_species_distances.c $(pkg-config --cflags --libs --static curupixa)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <curupixa/spr_distance.h>
#include <curupixa/reconciliation_roots.h>
#include <curupixa/taxon_map.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
  crpx_tree *gene=NULL,*species=NULL;
  crpx_ancestry_index *index=NULL;
  crpx_spr_workspace *workspace=NULL;
  crpx_root_reconciliations *roots=NULL;
  crpx_node_id *map=NULL;
  crpx_status status=crpx_tree_from_newick("(((A|1,A|2),C|1),(B|1,D|1));",&gene,NULL);
  if (status!=CRPX_OK) goto done;
  status=crpx_tree_from_newick("((A,B),(C,D));",&species,NULL);
  if (status!=CRPX_OK) goto done;
  {
    size_t nodes=crpx_tree_node_count(gene);
    map=(crpx_node_id *)malloc(nodes*sizeof(*map));
    if (!map) { status=CRPX_ERROR_ALLOCATION; goto done; }
    status=crpx_gene_species_map_labels(gene,species,'|',0,map,nodes);
    if (status!=CRPX_OK) goto done;
  }
  status=crpx_spr_workspace_new(&workspace);
  if (status==CRPX_OK) {
    crpx_spr_result distances;
    status=crpx_gene_species_spr_distance(gene,species,map,workspace,&distances);
    if (status==CRPX_OK)
      printf("RF=%zu dSPR=%zu crossing-prunes=%zu hdist=%" PRIu64 "\n",
             distances.rf,distances.spr,distances.spr_extra,distances.hdist);
  }
  if (status==CRPX_OK) status=crpx_ancestry_index_new(species,&index);
  if (status==CRPX_OK)
    status=crpx_reconcile_all_roots(gene,index,map,CRPX_RECONCILE_REPRESENTED_SPECIES,&roots);
  if (status==CRPX_OK) {
    size_t count,d,l,c;
    const crpx_root_reconciliation_score *scores=crpx_root_reconciliations_scores(roots,&count);
    status=crpx_root_reconciliations_minima(roots,&d,&l,&c);
    if (status==CRPX_OK)
      printf("%zu candidate edges: min duplications=%" PRIu64 ", losses=%" PRIu64 ", extra lineages=%" PRIu64 "\n",
             count,scores[d].duplications,scores[l].losses,scores[c].extra_lineages);
    /* d, l and c can identify different edges. Inspect scores to choose a
     * biologically appropriate joint criterion instead of combining minima. */
  }
done:
  if (status!=CRPX_OK) fprintf(stderr,"curupixa: %s\n",crpx_status_string(status));
  crpx_root_reconciliations_free(roots);
  crpx_ancestry_index_free(index);
  crpx_spr_workspace_free(workspace);
  free(map); crpx_tree_free(gene); crpx_tree_free(species);
  return status==CRPX_OK?EXIT_SUCCESS:EXIT_FAILURE;
}
