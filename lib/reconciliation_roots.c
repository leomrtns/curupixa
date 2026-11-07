/* Two-direction reconciliation derived from guenomu/genefam-dist topology_mrca.
 * Copyright (C) 2009-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "reconciliation_roots.h"
#include <stdlib.h>

typedef struct { crpx_node_id map; uint64_t duplications,losses,reduced; } message;
typedef struct { crpx_node_id neighbour[3]; uint8_t degree; message outgoing[3]; } vertex;
struct crpx_root_reconciliations { size_t count; crpx_root_reconciliation_score *scores; };
static bool add(uint64_t a,uint64_t b,uint64_t *out)
{ if (a>UINT64_MAX-b) return false; *out=a+b; return true; }
static crpx_status join(const crpx_ancestry_index *species,const uint32_t *reduced,
                        message a,message b,message *out)
{
  crpx_node_id p; crpx_status status=crpx_ancestry_index_lca(species,a.map,b.map,&p);
  if (status!=CRPX_OK) return status;
  bool duplication=p==a.map || p==b.map;
  uint64_t path=(uint64_t)(crpx_ancestry_index_depth(species,a.map)-crpx_ancestry_index_depth(species,p))+
    crpx_ancestry_index_depth(species,b.map)-crpx_ancestry_index_depth(species,p);
  uint64_t reduced_path=(uint64_t)(reduced[a.map]-reduced[p])+reduced[b.map]-reduced[p];
  if ((!duplication && (path<2 || reduced_path<2))) return CRPX_ERROR_TOPOLOGY;
  message r={p,0,0,0};
  if (!add(a.duplications,b.duplications,&r.duplications) ||
      !add(r.duplications,duplication,&r.duplications) ||
      !add(a.losses,b.losses,&r.losses) || !add(r.losses,path-(duplication?0:2),&r.losses) ||
      !add(a.reduced,b.reduced,&r.reduced) || !add(r.reduced,reduced_path-(duplication?0:2),&r.reduced))
    return CRPX_ERROR_OVERFLOW;
  *out=r; return CRPX_OK;
}
static size_t slot(const vertex *v,crpx_node_id node)
{ for (size_t i=0;i<v->degree;++i) if (v->neighbour[i]==node) return i; return SIZE_MAX; }
static crpx_status connect(vertex *v,crpx_node_id a,crpx_node_id b)
{
  if (v[a].degree>=3 || v[b].degree>=3) return CRPX_ERROR_UNSUPPORTED;
  v[a].neighbour[v[a].degree++]=b; v[b].neighbour[v[b].degree++]=a; return CRPX_OK;
}
static crpx_status component(vertex *v,crpx_node_id node,size_t excluded,const crpx_node_id *map,
                             const crpx_ancestry_index *species,const uint32_t *depth,message *out)
{
  if (v[node].degree==1) { *out=(message){map[node],0,0,0}; return CRPX_OK; }
  if (v[node].degree!=3) return CRPX_ERROR_TOPOLOGY;
  message inputs[2]; size_t count=0;
  for (size_t i=0;i<v[node].degree;++i) if (i!=excluded) {
    crpx_node_id other=v[node].neighbour[i]; size_t reverse=slot(v+other,node);
    if (reverse==SIZE_MAX) return CRPX_ERROR_TOPOLOGY;
    inputs[count++]=v[other].outgoing[reverse];
  }
  return join(species,depth,inputs[0],inputs[1],out);
}
crpx_status crpx_reconcile_all_roots(const crpx_tree *gene,const crpx_ancestry_index *species,
                                    const crpx_node_id *map,crpx_reconciliation_scope scope,
                                    crpx_root_reconciliations **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=NULL;
  if (!gene || !species || !map || (scope!=CRPX_RECONCILE_ALL_SPECIES &&
      scope!=CRPX_RECONCILE_REPRESENTED_SPECIES)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t gn=crpx_tree_node_count(gene),sn=crpx_ancestry_index_nodes(species);
  if (!gn || !crpx_ancestry_index_is_binary(species)) return CRPX_ERROR_UNSUPPORTED;
  if (gn>SIZE_MAX/sizeof(vertex) || gn>SIZE_MAX/sizeof(crpx_root_reconciliation_score) ||
      sn==SIZE_MAX || sn+1>SIZE_MAX/sizeof(size_t) || sn>SIZE_MAX/sizeof(uint64_t) ||
      gn>SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_status status=crpx_tree_validate(gene); if (status!=CRPX_OK) return status;
  uint64_t *presence=calloc(sn,sizeof(*presence)); uint32_t *depth=calloc(sn,sizeof(*depth));
  crpx_node_id *species_order=malloc(sn*sizeof(crpx_node_id));
  vertex *v=calloc(gn,sizeof(*v)); crpx_node_id *parent=malloc(gn*sizeof(*parent)),*order=malloc(gn*sizeof(*order));
  crpx_root_reconciliations *result=calloc(1,sizeof(*result));
  if (!presence || !depth || !species_order || !v || !parent || !order || !result) {
    status=CRPX_ERROR_ALLOCATION; goto done;
  }
  result->scores=calloc(gn,sizeof(*result->scores));
  if (!result->scores) { status=CRPX_ERROR_ALLOCATION; goto done; }
  size_t leaves=0,represented=0;
  for (size_t node=0;node<gn;++node) {
    crpx_node_id child=crpx_tree_first_child(gene,(crpx_node_id)node);
    if (child==CRPX_NODE_NONE) {
      if (!crpx_ancestry_index_is_leaf(species,map[node])) { status=CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      if (!presence[map[node]]) ++represented;
      ++presence[map[node]]; ++leaves;
    } else {
      crpx_node_id second=crpx_tree_next_sibling(gene,child);
      if (second==CRPX_NODE_NONE || crpx_tree_next_sibling(gene,second)!=CRPX_NODE_NONE) {
        status=CRPX_ERROR_UNSUPPORTED; goto done;
      }
    }
  }
  /* Counting-sort the species snapshot by depth, avoiding a borrowed traversal
   * or a quadratic scan. Parent depth is strictly smaller than child depth. */
  size_t *counts=calloc(sn+1,sizeof(*counts));
  if (!counts) { status=CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i=0;i<sn;++i) ++counts[crpx_ancestry_index_depth(species,(crpx_node_id)i)+1];
  for (size_t i=1;i<=sn;++i) counts[i]+=counts[i-1];
  for (size_t i=0;i<sn;++i) species_order[counts[crpx_ancestry_index_depth(species,(crpx_node_id)i)]++]=(crpx_node_id)i;
  free(counts);
  for (size_t i=sn;i>0;) {
    crpx_node_id node=species_order[--i],p=crpx_ancestry_index_parent(species,node);
    if (p!=CRPX_NODE_NONE) presence[p]+=presence[node];
  }
  /* Induced depth increases only at a split whose two descendant sides have
   * representatives, and at represented tips. */
  uint8_t *active_children=calloc(sn,sizeof(*active_children));
  if (!active_children) { status=CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t i=0;i<sn;++i) if (presence[i]) {
    crpx_node_id p=crpx_ancestry_index_parent(species,(crpx_node_id)i);
    if (p!=CRPX_NODE_NONE) ++active_children[p];
  }
  for (size_t i=0;i<sn;++i) {
    crpx_node_id node=species_order[i],p=crpx_ancestry_index_parent(species,node);
    depth[node]=(p==CRPX_NODE_NONE?0:depth[p])+
      (presence[node] && (crpx_ancestry_index_is_leaf(species,node) || active_children[node]==2));
  }
  free(active_children);
  if (leaves==1) { result->count=1; result->scores[0].first=result->scores[0].second=CRPX_NODE_NONE; goto success; }
  crpx_node_id root=crpx_tree_root(gene),left=crpx_tree_first_child(gene,root);
  crpx_node_id right=crpx_tree_next_sibling(gene,left);
  for (size_t node=0;node<gn;++node) {
    if (node==root) continue;
    crpx_node_id p=crpx_tree_parent(gene,(crpx_node_id)node);
    if (p!=root) { status=connect(v,(crpx_node_id)node,p); if (status!=CRPX_OK) goto done; }
  }
  status=connect(v,left,right); if (status!=CRPX_OK) goto done;
  /* Degree is at most three: keep candidate ordering independent of child order. */
  for (size_t node=0;node<gn;++node) for (size_t i=1;i<v[node].degree;++i) {
    crpx_node_id other=v[node].neighbour[i]; size_t j=i;
    while (j && v[node].neighbour[j-1]>other) { v[node].neighbour[j]=v[node].neighbour[j-1]; --j; }
    v[node].neighbour[j]=other;
  }
  for (size_t i=0;i<gn;++i) parent[i]=CRPX_NODE_NONE;
  size_t visited=1; order[0]=left;
  for (size_t i=0;i<visited;++i) {
    crpx_node_id node=order[i];
    for (size_t k=0;k<v[node].degree;++k) {
      crpx_node_id other=v[node].neighbour[k]; if (other==parent[node]) continue;
      parent[other]=node; order[visited++]=other;
    }
  }
  if (visited!=gn-1) { status=CRPX_ERROR_TOPOLOGY; goto done; }
  for (size_t i=visited;i>1;) {
    crpx_node_id node=order[--i]; size_t excluded=slot(v+node,parent[node]);
    status=component(v,node,excluded,map,species,depth,&v[node].outgoing[excluded]);
    if (status!=CRPX_OK) goto done;
  }
  for (size_t i=0;i<visited;++i) {
    crpx_node_id node=order[i];
    for (size_t k=0;k<v[node].degree;++k) if (v[node].neighbour[k]!=parent[node]) {
      status=component(v,node,k,map,species,depth,&v[node].outgoing[k]);
      if (status!=CRPX_OK) goto done;
    }
  }
  for (size_t node=0;node<gn;++node) for (size_t k=0;k<v[node].degree;++k) {
    crpx_node_id other=v[node].neighbour[k]; if (node>=other) continue;
    message r; status=join(species,depth,v[node].outgoing[k],v[other].outgoing[slot(v+other,(crpx_node_id)node)],&r);
    if (status!=CRPX_OK) goto done;
    uint64_t difference=2*(uint64_t)(leaves-represented),sum;
    if (!add(r.reduced,difference,&sum) || r.duplications>UINT64_MAX/2 || sum<2*r.duplications) {
      status=CRPX_ERROR_OVERFLOW; goto done;
    }
    result->scores[result->count++]=(crpx_root_reconciliation_score){(crpx_node_id)node,other,r.duplications,
      scope==CRPX_RECONCILE_ALL_SPECIES?r.losses:r.reduced,sum-2*r.duplications};
  }
success:
  *out=result; result=NULL;
done:
  crpx_root_reconciliations_free(result); free(presence); free(depth); free(species_order); free(v); free(parent); free(order);
  return status;
}
void crpx_root_reconciliations_free(crpx_root_reconciliations *r) { if (r) { free(r->scores); free(r); } }
const crpx_root_reconciliation_score *crpx_root_reconciliations_scores(const crpx_root_reconciliations *r,size_t *count)
{ if (count) *count=r?r->count:0; return r?r->scores:NULL; }
crpx_status crpx_root_reconciliations_minima(const crpx_root_reconciliations *r,size_t *dups,size_t *loss,size_t *coal)
{
  if (!r || !dups || !loss || !coal || !r->count) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t d=0,l=0,c=0;
  for (size_t i=1;i<r->count;++i) {
    if (r->scores[i].duplications<r->scores[d].duplications) d=i;
    if (r->scores[i].losses<r->scores[l].losses) l=i;
    if (r->scores[i].extra_lineages<r->scores[c].extra_lineages) c=i;
  }
  *dups=d; *loss=l; *coal=c; return CRPX_OK;
}
