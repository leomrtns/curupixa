/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "ancestry_index.h"
#include <stdlib.h>
struct crpx_ancestry_index {
  size_t nodes,levels;
  crpx_node_id *ancestors;
  uint32_t *depth;
  uint8_t *leaves;
  bool binary;
};
crpx_status crpx_ancestry_index_new(const crpx_tree *tree,crpx_ancestry_index **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  size_t n=crpx_tree_node_count(tree);
  if (!tree || !n) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_tree_workspace *w=NULL;
  crpx_status status=crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL|CRPX_TREE_WORKSPACE_DEPTH,&w);
  if (status==CRPX_OK) status=crpx_tree_workspace_refresh(w,tree);
  if (status!=CRPX_OK) { crpx_tree_workspace_free(w); return status; }
  size_t levels=1; for (size_t x=n;x>1;x>>=1) ++levels;
  if (n>SIZE_MAX/levels/sizeof(crpx_node_id) || n>SIZE_MAX/sizeof(uint32_t)) {
    crpx_tree_workspace_free(w); return CRPX_ERROR_OVERFLOW;
  }
  crpx_ancestry_index *index=calloc(1,sizeof(*index));
  if (!index) { crpx_tree_workspace_free(w); return CRPX_ERROR_ALLOCATION; }
  index->nodes=n; index->levels=levels; index->binary=true;
  index->ancestors=malloc(n*levels*sizeof(crpx_node_id));
  index->depth=malloc(n*sizeof(uint32_t)); index->leaves=calloc(n,sizeof(uint8_t));
  if (!index->ancestors || !index->depth || !index->leaves) {
    crpx_tree_workspace_free(w); crpx_ancestry_index_free(index); return CRPX_ERROR_ALLOCATION;
  }
  const uint32_t *depth=crpx_tree_workspace_depth(w);
  for (size_t i=0;i<n;++i) {
    index->depth[i]=depth[i]; index->ancestors[i]=crpx_tree_parent(tree,(crpx_node_id)i);
    size_t children=0;
    for (crpx_node_id child=crpx_tree_first_child(tree,(crpx_node_id)i);child!=CRPX_NODE_NONE;
         child=crpx_tree_next_sibling(tree,child)) ++children;
    index->leaves[i]=children==0;
    if (children && children!=2) index->binary=false;
  }
  for (size_t level=1;level<levels;++level) for (size_t i=0;i<n;++i) {
    crpx_node_id parent=index->ancestors[(level-1)*n+i];
    index->ancestors[level*n+i]=parent==CRPX_NODE_NONE?CRPX_NODE_NONE:index->ancestors[(level-1)*n+parent];
  }
  crpx_tree_workspace_free(w); *out=index; return CRPX_OK;
}
void crpx_ancestry_index_free(crpx_ancestry_index *i)
{ if (i) { free(i->ancestors); free(i->depth); free(i->leaves); free(i); } }
size_t crpx_ancestry_index_nodes(const crpx_ancestry_index *i) { return i?i->nodes:0; }
size_t crpx_ancestry_index_storage_size(const crpx_ancestry_index *i)
{ return i?sizeof(*i)+i->nodes*(i->levels*sizeof(crpx_node_id)+sizeof(uint32_t)+sizeof(uint8_t)):0; }
crpx_node_id crpx_ancestry_index_parent(const crpx_ancestry_index *i,crpx_node_id n)
{ return i && n<i->nodes?i->ancestors[n]:CRPX_NODE_NONE; }
uint32_t crpx_ancestry_index_depth(const crpx_ancestry_index *i,crpx_node_id n)
{ return i && n<i->nodes?i->depth[n]:UINT32_MAX; }
bool crpx_ancestry_index_is_leaf(const crpx_ancestry_index *i,crpx_node_id n)
{ return i && n<i->nodes && i->leaves[n]; }
bool crpx_ancestry_index_is_binary(const crpx_ancestry_index *i) { return i && i->binary; }
crpx_status crpx_ancestry_index_lca(const crpx_ancestry_index *i,crpx_node_id a,crpx_node_id b,crpx_node_id *out)
{
  if (!i || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (a>=i->nodes || b>=i->nodes) return CRPX_ERROR_OUT_OF_RANGE;
  if (i->depth[a]<i->depth[b]) { crpx_node_id swap=a; a=b; b=swap; }
  uint32_t difference=i->depth[a]-i->depth[b];
  for (size_t level=0;difference;++level,difference>>=1) if (difference&1) a=i->ancestors[level*i->nodes+a];
  if (a==b) { *out=a; return CRPX_OK; }
  for (size_t level=i->levels;level>0;) {
    --level;
    crpx_node_id x=i->ancestors[level*i->nodes+a],y=i->ancestors[level*i->nodes+b];
    if (x!=y) { a=x; b=y; }
  }
  *out=i->ancestors[a]; return CRPX_OK;
}
