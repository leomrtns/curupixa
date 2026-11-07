/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "taxon_map.h"
#include "hashtable.h"
#include <stdlib.h>
#include <string.h>
crpx_status crpx_gene_species_map_labels(const crpx_tree *gene,const crpx_tree *species,
                                        char delimiter,size_t field,crpx_node_id *out,size_t capacity)
{
  if (!gene || !species || !out || (!delimiter && field)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t n=crpx_tree_node_count(gene);
  if (!n || !crpx_tree_node_count(species)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity<n) return CRPX_ERROR_OUT_OF_RANGE;
  if (n>SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *map=malloc(n*sizeof(*map)); if (!map) return CRPX_ERROR_ALLOCATION;
  crpx_hashtable_t labels=NULL; crpx_status status=crpx_hashtable_new(crpx_tree_leaf_count(species),0,&labels);
  for (size_t i=0;i<crpx_tree_node_count(species) && status==CRPX_OK;++i)
    if (crpx_tree_first_child(species,(crpx_node_id)i)==CRPX_NODE_NONE) {
      const char *label=crpx_tree_label(species,(crpx_node_id)i); bool replaced;
      if (!label || !*label) { status=CRPX_ERROR_INVALID_ARGUMENT; break; }
      status=crpx_hashtable_put_string(labels,label,i,&replaced);
      if (status==CRPX_OK && replaced) status=CRPX_ERROR_INVALID_ARGUMENT;
    }
  for (size_t i=0;i<n && status==CRPX_OK;++i) {
    map[i]=CRPX_NODE_NONE;
    if (crpx_tree_first_child(gene,(crpx_node_id)i)!=CRPX_NODE_NONE) continue;
    const char *name=crpx_tree_label(gene,(crpx_node_id)i);
    if (!name) { status=CRPX_ERROR_INVALID_ARGUMENT; break; }
    for (size_t f=0;f<field;++f) {
      name=strchr(name,delimiter); if (!name) break; ++name;
    }
    if (!name || !*name) { status=CRPX_ERROR_NOT_FOUND; break; }
    const char *end=delimiter?strchr(name,delimiter):NULL;
    size_t length=end?(size_t)(end-name):strlen(name); uint64_t id;
    if (!length) { status=CRPX_ERROR_NOT_FOUND; break; }
    status=crpx_hashtable_get(labels,name,length,&id);
    if (status==CRPX_OK) map[i]=(crpx_node_id)id;
  }
  if (status==CRPX_OK) memcpy(out,map,n*sizeof(*map));
  free(map); crpx_hashtable_free(labels); return status;
}
