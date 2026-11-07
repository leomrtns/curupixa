/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "tree.h"
#include "topology_distance.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OK(x) assert((x)==CRPX_OK)
static void edge_cases(const char *text, crpx_tree_edge_root_policy policy)
{
  crpx_tree *source; OK(crpx_tree_from_newick(text,&source,NULL));
  size_t n=crpx_tree_node_count(source);
  for (size_t child=0;child<n;++child) if (child!=crpx_tree_root(source)) {
    for (size_t f=0;f<3;++f) {
      crpx_tree *tree; OK(crpx_tree_clone(source,&tree));
      crpx_tree_workspace *workspace; OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL,&workspace));
      OK(crpx_tree_workspace_refresh(workspace,tree));
      crpx_node_id root; uint64_t revision=crpx_tree_revision(tree);
      OK(crpx_tree_reroot_edge(tree,(crpx_node_id)child,0.5*f,policy,&root));
      assert(crpx_tree_revision(tree)==revision+1 && !crpx_tree_workspace_is_current(workspace,tree));
      assert(root==crpx_tree_root(tree));
      assert(crpx_tree_node_count(tree)==n+(policy==CRPX_TREE_KEEP_OLD_ROOT));
      assert(crpx_tree_leaf_count(tree)==crpx_tree_leaf_count(source));
      OK(crpx_tree_validate(tree));
      double stem; OK(crpx_tree_branch_length(tree,root,&stem)); assert(stem==23);
      for (size_t i=0;i<n;++i) {
        const char *a=crpx_tree_label(source,(crpx_node_id)i),*b=crpx_tree_label(tree,(crpx_node_id)i);
        assert((!a && !b) || (a && b && !strcmp(a,b)));
        if (crpx_tree_first_child(source,(crpx_node_id)i)!=CRPX_NODE_NONE) continue;
        assert(crpx_tree_first_child(tree,(crpx_node_id)i)==CRPX_NODE_NONE);
        for (size_t j=0;j<n;++j) if (crpx_tree_first_child(source,(crpx_node_id)j)==CRPX_NODE_NONE) {
          double before,after;
          OK(crpx_tree_patristic_distance(source,(crpx_node_id)i,(crpx_node_id)j,NAN,&before));
          OK(crpx_tree_patristic_distance(tree,(crpx_node_id)i,(crpx_node_id)j,NAN,&after));
          assert(fabs(before-after)<1e-12);
        }
      }
      crpx_tree_workspace_free(workspace); crpx_tree_free(tree);
    }
  }
  crpx_tree_free(source);
}
static void root_policies(void)
{
  edge_cases("((A:1,B:2)X:3,(C:4,D:0)Y:5)R:23;",CRPX_TREE_REUSE_BINARY_ROOT);
  edge_cases("(A:1,(B:2,(C:3,D:4):5):6)R:23;",CRPX_TREE_REUSE_BINARY_ROOT);
  edge_cases("(A:1,B:2)R:23;",CRPX_TREE_REUSE_BINARY_ROOT);
  edge_cases("((A:1,B:2):3,C:4,D:0)R:23;",CRPX_TREE_KEEP_OLD_ROOT);
  edge_cases("((A:1,B:2):3,(C:4,D:0):5)R:23;",CRPX_TREE_KEEP_OLD_ROOT);
  edge_cases("(A:-1,(B:2,C:3):-4)R:23;",CRPX_TREE_REUSE_BINARY_ROOT);
  crpx_tree *tree; OK(crpx_tree_from_newick("(A:1,B:2,C:7)R:23;",&tree,NULL));
  crpx_node_id result=77; char *before,*after; OK(crpx_tree_to_newick(tree,&before));
  assert(crpx_tree_reroot_edge(tree,3,.5,CRPX_TREE_REUSE_BINARY_ROOT,&result)==CRPX_ERROR_UNSUPPORTED);
  assert(result==77); OK(crpx_tree_to_newick(tree,&after)); assert(!strcmp(before,after));
  free(before); free(after);
  /* B--C diameter is 9; its midpoint is 4.5 from C, or 2.5 from C's parent. */
  OK(crpx_tree_reroot_edge(tree,3,2.5/7,CRPX_TREE_KEEP_OLD_ROOT,&result));
  double a,b; OK(crpx_tree_patristic_distance(tree,result,2,NAN,&a));
  OK(crpx_tree_patristic_distance(tree,result,3,NAN,&b)); assert(fabs(a-4.5)<1e-12 && fabs(b-4.5)<1e-12);
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(((A:1,B:2):3):4):5;",&tree,NULL));
  assert(crpx_tree_reroot_edge(tree,3,.5,CRPX_TREE_KEEP_OLD_ROOT,&result)==CRPX_ERROR_UNSUPPORTED);
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A,B);",&tree,NULL));
  OK(crpx_tree_reroot_edge(tree,1,.5,CRPX_TREE_REUSE_BINARY_ROOT,&result));
  assert(crpx_tree_branch_length(tree,1,&a)==CRPX_ERROR_NOT_FOUND); OK(crpx_tree_validate(tree));
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A:1e308,B:1e308);",&tree,NULL));
  assert(crpx_tree_reroot_edge(tree,1,.5,CRPX_TREE_REUSE_BINARY_ROOT,&result)==CRPX_ERROR_OVERFLOW);
  crpx_tree_free(tree);
}
static void build_stars(void)
{
  crpx_tree_builder *builder; OK(crpx_tree_builder_new(&builder));
  for (size_t size=1000;size<=4000;size*=2) {
    crpx_tree *tree; OK(crpx_tree_new(&tree)); OK(crpx_tree_reserve(tree,size+2));
    crpx_node_id root,node; OK(crpx_tree_builder_append(builder,tree,CRPX_NODE_NONE,&root));
    for (size_t i=0;i<size;++i) {
      OK(crpx_tree_builder_append(builder,tree,root,&node)); assert(node==i+1);
      OK(crpx_tree_set_label(tree,node,"tip")); OK(crpx_tree_set_branch_length(tree,node,1));
    }
    assert(crpx_tree_leaf_count(tree)==size);
    for (size_t i=1;i<size;++i) assert(crpx_tree_next_sibling(tree,(crpx_node_id)i)==i+1);
    OK(crpx_tree_set_label(tree,root,"root"));
    OK(crpx_tree_builder_append(builder,tree,root,&node)); assert(node==size+1);
    OK(crpx_tree_validate(tree));
    /* Tail allocation belongs to the builder; the tree retains exactly 12 bytes/node. */
    assert(crpx_tree_builder_storage_size(builder)>size*sizeof(crpx_node_id));
    crpx_tree_builder_reset(builder); crpx_tree_free(tree);
  }
  crpx_tree_builder_free(builder);
}
static void validation(void)
{
  const uint32_t names=CRPX_TREE_REQUIRE_NAMED_TIPS|CRPX_TREE_REQUIRE_UNIQUE_TIPS;
  const uint32_t lengths=CRPX_TREE_REQUIRE_COMPLETE_LENGTHS|CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS;
  crpx_tree *tree; OK(crpx_tree_from_newick("((A:0,B:2):3,C:4):-10;",&tree,NULL));
  OK(crpx_tree_validate_policy(tree,names|lengths|CRPX_TREE_REQUIRE_BINARY));
  assert(crpx_tree_validate_policy(tree,lengths|CRPX_TREE_INCLUDE_ROOT_STEM)==CRPX_ERROR_INVALID_ARGUMENT);
  assert(crpx_tree_validate_policy(tree,CRPX_TREE_REQUIRE_POSITIVE_LENGTHS)==CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A,A,);",&tree,NULL));
  OK(crpx_tree_validate_policy(tree,0));
  assert(crpx_tree_validate_policy(tree,CRPX_TREE_REQUIRE_NAMED_TIPS)==CRPX_ERROR_INVALID_ARGUMENT);
  assert(crpx_tree_validate_policy(tree,CRPX_TREE_REQUIRE_UNIQUE_TIPS)==CRPX_ERROR_INVALID_ARGUMENT);
  assert(crpx_tree_validate_policy(tree,CRPX_TREE_REQUIRE_COMPLETE_LENGTHS)==CRPX_ERROR_NOT_FOUND);
  assert(crpx_tree_validate_policy(tree,CRPX_TREE_REQUIRE_BINARY)==CRPX_ERROR_UNSUPPORTED);
  assert(crpx_tree_validate_policy(tree,UINT32_C(128))==CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
}
int main(void)
{
  root_policies(); build_stars(); validation();
  puts("edge rooting, stable IDs/distances, optional linear builder and policy validation: passed");
  return 0;
}
