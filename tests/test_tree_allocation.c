/* Fault injection compiles the implementation locally so allocator interception
 * remains confined to this test; production code needs no allocator hooks.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
static size_t allocation_calls, fail_call;
static void *checked_malloc(size_t size)
{ return ++allocation_calls==fail_call ? NULL : malloc(size); }
static void *checked_calloc(size_t count,size_t size)
{ return ++allocation_calls==fail_call ? NULL : calloc(count,size); }
static void *checked_realloc(void *ptr,size_t size)
{ return ++allocation_calls==fail_call ? NULL : realloc(ptr,size); }
#define malloc checked_malloc
#define calloc checked_calloc
#define realloc checked_realloc
#include "../lib/tree.c"
#undef malloc
#undef calloc
#undef realloc
#define OK(x) assert((x)==CRPX_OK)
int main(void)
{
  for (unsigned policy=0;policy<2;++policy) {
    size_t failures=0;
    for (size_t fail=1;fail<16;++fail) {
      crpx_tree *tree; OK(crpx_tree_from_newick("((A:1,B:2):3,(C:4,D:5):6)R:7;",&tree,NULL));
      char *before,*after; OK(crpx_tree_to_newick(tree,&before));
      uint64_t revision=crpx_tree_revision(tree); size_t bytes=crpx_tree_storage_size(tree);
      allocation_calls=0; fail_call=fail; crpx_node_id root=91;
      crpx_status status=crpx_tree_reroot_edge(tree,3,.4,(crpx_tree_edge_root_policy)policy,&root);
      fail_call=0;
      if (status==CRPX_ERROR_ALLOCATION) {
        ++failures; assert(root==91 && revision==crpx_tree_revision(tree));
        assert(bytes==crpx_tree_storage_size(tree));
        OK(crpx_tree_to_newick(tree,&after)); assert(!strcmp(before,after)); free(after);
      } else { OK(status); free(before); crpx_tree_free(tree); break; }
      free(before); crpx_tree_free(tree);
    }
    assert(failures==5); /* validate:2, replacement connectivity/lengths/labels:3 */
  }
  for (size_t fail=1;fail<=2;++fail) {
    crpx_tree *tree; crpx_tree_builder *builder; crpx_node_id node=93;
    OK(crpx_tree_new(&tree)); OK(crpx_tree_builder_new(&builder));
    allocation_calls=0; fail_call=fail;
    assert(crpx_tree_builder_append(builder,tree,CRPX_NODE_NONE,&node)==CRPX_ERROR_ALLOCATION);
    fail_call=0; assert(node==93 && !crpx_tree_node_count(tree));
    OK(crpx_tree_builder_append(builder,tree,CRPX_NODE_NONE,&node)); assert(node==0);
    crpx_tree_builder_free(builder); crpx_tree_free(tree);
  }
  /* Prove annotation between appends does not trigger O(V) refresh: with spare
   * capacities, each checked append performs zero allocations, whereas a refresh
   * would allocate the structural validator's arrays. */
  for (size_t size=1000;size<=4000;size*=2) {
    crpx_tree *tree; crpx_tree_builder *builder; crpx_node_id root,node=CRPX_NODE_NONE;
    OK(crpx_tree_new(&tree)); OK(crpx_tree_reserve(tree,size+1));
    OK(crpx_tree_builder_new(&builder));
    OK(crpx_tree_builder_append(builder,tree,CRPX_NODE_NONE,&root));
    for (size_t i=0;i<size;++i) {
      /* Tail capacity itself may grow geometrically, but never from a refresh. */
      size_t old_capacity=builder->capacity;
      allocation_calls=0; OK(crpx_tree_builder_append(builder,tree,root,&node));
      assert(allocation_calls==(old_capacity==builder->capacity ? 0u : 1u));
      uint64_t topology=tree->topology_revision;
      OK(crpx_tree_set_label(tree,node,"tip")); OK(crpx_tree_set_branch_length(tree,node,1));
      assert(tree->topology_revision==topology && builder->revision==topology);
    }
    /* A genuine external topology edit does invalidate the tail snapshot. */
    OK(crpx_tree_reparent(tree,node,1));
    allocation_calls=0; OK(crpx_tree_builder_append(builder,tree,root,&node));
    assert(allocation_calls>=2); /* validator allocations prove refresh */
    crpx_tree_builder_free(builder); crpx_tree_free(tree);
  }
  puts("tree root-insertion and builder allocation failures remain recoverable: passed");
  return 0;
}
