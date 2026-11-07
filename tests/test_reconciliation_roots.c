/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "reconciliation_roots.h"
#include "reconciliation.h"
#include "taxon_map.h"
#include "topology_distance.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(x) assert((x)==CRPX_OK)
static crpx_node_id leaf(const crpx_tree *t,const char *name)
{
  for(size_t i=0;i<crpx_tree_node_count(t);++i) {
    const char *s=crpx_tree_label(t,(crpx_node_id)i);
    if(s && !strcmp(s,name) && crpx_tree_first_child(t,(crpx_node_id)i)==CRPX_NODE_NONE) return (crpx_node_id)i;
  }
  return CRPX_NODE_NONE;
}
static void test_index(void)
{
  crpx_tree *t; OK(crpx_tree_from_newick("(((A,B),C),(D,E));",&t,NULL));
  crpx_ancestry_index *idx; OK(crpx_ancestry_index_new(t,&idx));
  for(size_t a=0;a<crpx_tree_node_count(t);++a) for(size_t b=0;b<crpx_tree_node_count(t);++b) {
    crpx_node_id x,y; OK(crpx_tree_lca(t,(crpx_node_id)a,(crpx_node_id)b,&x));
    OK(crpx_ancestry_index_lca(idx,(crpx_node_id)a,(crpx_node_id)b,&y)); assert(x==y);
  }
  crpx_node_id a=leaf(t,"A"),b=leaf(t,"B"),expected=crpx_tree_parent(t,a),result;
  crpx_tree_free(t); OK(crpx_ancestry_index_lca(idx,a,b,&result)); assert(result==expected);
  crpx_ancestry_index_free(idx);
}
/* Test-only independent re-rooting rebuilds a fresh tree for each edge. */
typedef struct { crpx_node_id adjacent[32][3]; size_t degree[32]; } graph;
static void link_nodes(graph *g,crpx_node_id a,crpx_node_id b)
{ g->adjacent[a][g->degree[a]++]=b; g->adjacent[b][g->degree[b]++]=a; }
static void copy_component(const graph *g,const crpx_tree *source,crpx_tree *target,
                            crpx_node_id node,crpx_node_id previous,crpx_node_id parent)
{
  crpx_node_id added; OK(crpx_tree_add_node(target,parent,&added));
  if(crpx_tree_first_child(source,node)==CRPX_NODE_NONE) OK(crpx_tree_set_label(target,added,crpx_tree_label(source,node)));
  for(size_t i=0;i<g->degree[node];++i) if(g->adjacent[node][i]!=previous)
    copy_component(g,source,target,g->adjacent[node][i],node,added);
}
static char *induced(const crpx_tree *species,crpx_node_id node,const bool *present)
{
  crpx_node_id a=crpx_tree_first_child(species,node);
  if(a==CRPX_NODE_NONE) {
    if(!present[node]) return NULL;
    const char *s=crpx_tree_label(species,node); char *text=malloc(strlen(s)+1); assert(text); strcpy(text,s); return text;
  }
  crpx_node_id b=crpx_tree_next_sibling(species,a);
  char *x=induced(species,a,present),*y=induced(species,b,present);
  if(!x) return y;
  if(!y) return x;
  char *r=malloc(strlen(x)+strlen(y)+4); assert(r); sprintf(r,"(%s,%s)",x,y); free(x); free(y); return r;
}
static void case_roots(const char *newick,const char *species_newick)
{
  crpx_tree *gene,*species,*reduced;
  OK(crpx_tree_from_newick(newick,&gene,NULL));
  OK(crpx_tree_from_newick(species_newick,&species,NULL));
  size_t n=crpx_tree_node_count(gene); assert(n<32);
  crpx_node_id map[32]; OK(crpx_gene_species_map_labels(gene,species,'|',0,map,32));
  bool present[32]={0}; assert(crpx_tree_node_count(species)<32);
  for(size_t i=0;i<n;++i) if(map[i]!=CRPX_NODE_NONE) present[map[i]]=true;
  char *text=induced(species,crpx_tree_root(species),present); assert(text);
  char *terminated=malloc(strlen(text)+2); assert(terminated); sprintf(terminated,"%s;",text);
  OK(crpx_tree_from_newick(terminated,&reduced,NULL)); free(text); free(terminated);
  crpx_ancestry_index *index; OK(crpx_ancestry_index_new(species,&index));
  crpx_node_id root=crpx_tree_root(gene),left=crpx_tree_first_child(gene,root);
  graph g={0};
  for(size_t i=0;i<n;++i) if(i!=root) {
    crpx_node_id p=crpx_tree_parent(gene,(crpx_node_id)i); if(p!=root) link_nodes(&g,(crpx_node_id)i,p);
  }
  link_nodes(&g,left,crpx_tree_next_sibling(gene,left));
  for(unsigned scope=0;scope<2;++scope) {
    crpx_root_reconciliations *r; OK(crpx_reconcile_all_roots(gene,index,map,(crpx_reconciliation_scope)scope,&r));
    size_t count; const crpx_root_reconciliation_score *scores=crpx_root_reconciliations_scores(r,&count);
    assert(count==2*crpx_tree_leaf_count(gene)-3);
    for(size_t i=0;i<count;++i) {
      crpx_tree *oriented; crpx_node_id newroot;
      OK(crpx_tree_new(&oriented)); OK(crpx_tree_add_node(oriented,CRPX_NODE_NONE,&newroot));
      copy_component(&g,gene,oriented,scores[i].first,scores[i].second,newroot);
      copy_component(&g,gene,oriented,scores[i].second,scores[i].first,newroot);
      crpx_node_id mapped[32]; crpx_reconciliation *expected,*coal;
      const crpx_tree *loss_species=scope?reduced:species;
      OK(crpx_gene_species_map_labels(oriented,loss_species,'|',0,mapped,32));
      OK(crpx_reconcile_lca(oriented,loss_species,mapped,&expected));
      assert(scores[i].duplications==crpx_reconciliation_duplications(expected));
      assert(scores[i].losses==crpx_reconciliation_losses(expected));
      OK(crpx_gene_species_map_labels(oriented,reduced,'|',0,mapped,32));
      OK(crpx_reconcile_lca(oriented,reduced,mapped,&coal));
      assert(scores[i].extra_lineages==crpx_reconciliation_extra_lineages(coal));
      crpx_reconciliation_free(expected); crpx_reconciliation_free(coal); crpx_tree_free(oriented);
    }
    size_t d,l,c; OK(crpx_root_reconciliations_minima(r,&d,&l,&c));
    for(size_t i=0;i<count;++i) assert(scores[d].duplications<=scores[i].duplications &&
      scores[l].losses<=scores[i].losses && scores[c].extra_lineages<=scores[i].extra_lineages);
    crpx_root_reconciliations_free(r);
  }
  crpx_node_id bad[32]; memcpy(bad,map,n*sizeof(*bad)); bad[leaf(gene,"A|1")]=CRPX_NODE_NONE;
  crpx_root_reconciliations *r=(void *)1;
  assert(crpx_reconcile_all_roots(gene,index,bad,CRPX_RECONCILE_ALL_SPECIES,&r)==CRPX_ERROR_INVALID_ARGUMENT && !r);
  crpx_ancestry_index_free(index); crpx_tree_free(gene); crpx_tree_free(species); crpx_tree_free(reduced);
}
static void test_map(void)
{
  crpx_tree *g,*s; OK(crpx_tree_from_newick("(copy|A,copy|B);",&g,NULL));
  OK(crpx_tree_from_newick("(A,B);",&s,NULL)); crpx_node_id map[3];
  OK(crpx_gene_species_map_labels(g,s,'|',1,map,3)); assert(map[1]==leaf(s,"A"));
  crpx_node_id before[3]; memcpy(before,map,sizeof(map));
  assert(crpx_gene_species_map_labels(g,s,'|',2,map,3)==CRPX_ERROR_NOT_FOUND);
  assert(!memcmp(before,map,sizeof(map))); crpx_tree_free(g); crpx_tree_free(s);
}
static void test_single_and_polytomy(void)
{
  crpx_tree *g,*s; crpx_ancestry_index *index; crpx_root_reconciliations *r;
  OK(crpx_tree_from_newick("A;",&g,NULL)); OK(crpx_tree_from_newick("((A,B),C);",&s,NULL));
  OK(crpx_ancestry_index_new(s,&index)); crpx_node_id map[1];
  OK(crpx_gene_species_map_labels(g,s,0,0,map,1));
  OK(crpx_reconcile_all_roots(g,index,map,CRPX_RECONCILE_ALL_SPECIES,&r));
  crpx_tree_free(g); crpx_ancestry_index_free(index); crpx_tree_free(s);
  size_t n; const crpx_root_reconciliation_score *scores=crpx_root_reconciliations_scores(r,&n);
  assert(n==1 && scores[0].first==CRPX_NODE_NONE && scores[0].second==CRPX_NODE_NONE);
  assert(!scores[0].duplications && !scores[0].losses && !scores[0].extra_lineages);
  crpx_root_reconciliations_free(r);
  OK(crpx_tree_from_newick("(A,B,C);",&g,NULL)); OK(crpx_tree_from_newick("((A,B),C);",&s,NULL));
  OK(crpx_ancestry_index_new(s,&index)); crpx_node_id many[4];
  OK(crpx_gene_species_map_labels(g,s,0,0,many,4));
  assert(crpx_reconcile_all_roots(g,index,many,CRPX_RECONCILE_ALL_SPECIES,&r)==CRPX_ERROR_UNSUPPORTED && !r);
  crpx_tree_free(g); crpx_tree_free(s); crpx_ancestry_index_free(index);
}
int main(void)
{
  test_index(); test_map(); test_single_and_polytomy();
  const char *trees[]={"((A|1,B|1),(C|1,D|1));","((A|1,C|1),(B|1,D|1));",
    "(((A|1,A|2),C|1),(B|1,D|1));","((A|1,C|1),(A|2,C|2));","(A|1,C|1);",
    "(((A|1,A|2),A|3),A|4);"};
  for(size_t i=0;i<sizeof(trees)/sizeof(*trees);++i) {
    case_roots(trees[i],"((A,B),(C,D));");
    case_roots(trees[i],"((((A,B),E),C),D);"); /* absent side branches at unequal depths */
  }
  puts("indexed LCA and all-root reconciliation vs independent rooted reconstructions: passed"); return 0;
}
