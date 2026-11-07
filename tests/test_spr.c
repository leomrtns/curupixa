/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "spr_distance.h"
#include "taxon_map.h"
#include "topology_distance.h"
#include "topology_randomise.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(x) assert((x)==CRPX_OK)
static uint64_t rng(void *state)
{ uint64_t *x=state; *x=*x*UINT64_C(6364136223846793005)+1; return *x; }
static void pair(crpx_spr_workspace *w,const char *first,const char *second,size_t spr,size_t rf)
{
  crpx_tree *a,*b; OK(crpx_tree_from_newick(first,&a,NULL)); OK(crpx_tree_from_newick(second,&b,NULL));
  crpx_spr_result r; uint64_t revision=crpx_tree_revision(a);
  OK(crpx_tree_spr_distance(a,b,w,&r)); assert(r.spr==spr && r.rf==rf);
  assert(revision==crpx_tree_revision(a)); size_t exact; OK(crpx_tree_rf_distance(a,b,false,&exact)); assert(r.rf==exact);
  crpx_tree_free(a); crpx_tree_free(b);
}
int main(void)
{
  crpx_spr_workspace *w; OK(crpx_spr_workspace_new(&w)); size_t empty=crpx_spr_workspace_storage_size(w);
  pair(w,"((A,B),(C,D));","((A,B),(C,D));",0,0);
  pair(w,"((A,B),(C,D));","((A,C),(B,D));",1,2);
  pair(w,"((A,B),(C,D));","(A,B,C,D);",0,1); /* unresolved compatibility */
  pair(w,"(((A,B),C),(D,(E,F)));","(((A,C),B),(D,(E,F)));",1,2);
  uint64_t state=19;
  for(size_t size=4;size<=130;size=size<10?size+1:130) {
    crpx_tree *a,*b; OK(crpx_tree_random_binary(size,NULL,rng,&state,&a)); OK(crpx_tree_clone(a,&b));
    crpx_spr_result result; OK(crpx_tree_spr_distance(a,a,w,&result)); assert(!result.rf && !result.spr);
    for(size_t iteration=0;iteration<8;++iteration) {
      crpx_status move=crpx_tree_random_spr(b,.5,rng,&state);
      assert(move==CRPX_OK || move==CRPX_ERROR_NOT_FOUND);
      OK(crpx_tree_spr_distance(a,b,w,&result)); size_t rf;
      OK(crpx_tree_rf_distance(a,b,false,&rf)); assert(result.rf==rf && result.spr<=size-3);
      /* Both binary trees have equal split counts: positive/zero dummy policies
       * coincide, providing an independent h-distance implementation. */
      size_t nodes=crpx_tree_node_count(a); crpx_node_id *map=malloc(nodes*sizeof(*map)); assert(map);
      OK(crpx_gene_species_map_labels(a,b,0,0,map,nodes)); uint64_t hdist;
      OK(crpx_gene_species_split_distances(a,b,map,&rf,&hdist));
      assert(result.rf==rf && result.hdist==hdist); free(map);
    }
    crpx_tree_free(a); crpx_tree_free(b); if(size==130) break;
  }
  crpx_tree *gene,*species;
  OK(crpx_tree_from_newick("(((A|1,A|2),B|1),(C|1,D|1));",&gene,NULL));
  OK(crpx_tree_from_newick("((A,B),(C,D));",&species,NULL));
  crpx_node_id map[9]; OK(crpx_gene_species_map_labels(gene,species,'|',0,map,9));
  crpx_spr_result r; OK(crpx_gene_species_spr_distance(gene,species,map,w,&r)); assert(!r.rf && !r.spr);
  crpx_tree_free(gene); crpx_tree_free(species);
  OK(crpx_tree_from_newick("((A,B),(C,D));",&gene,NULL));
  OK(crpx_tree_from_newick("((A,C),(B,D));",&species,NULL));
  OK(crpx_tree_spr_distance(gene,species,w,&r));
  assert(r.hdist==2 && r.reduced_hdist==2); /* legacy skipped the one-pair assignment */
  crpx_tree_free(species); OK(crpx_tree_from_newick("((A,C),(B,E));",&species,NULL));
  crpx_spr_result before=r;
  assert(crpx_tree_spr_distance(gene,species,w,&r)==CRPX_ERROR_INVALID_ARGUMENT);
  assert(!memcmp(&r,&before,sizeof(r))); crpx_tree_free(gene); crpx_tree_free(species);
  assert(crpx_spr_workspace_storage_size(w)>empty); crpx_spr_workspace_clear(w);
  assert(crpx_spr_workspace_storage_size(w)==empty); crpx_spr_workspace_free(w);
  puts("dSPR split pruning, mapped copies, word boundaries and reusable workspace: passed"); return 0;
}
