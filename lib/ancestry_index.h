/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_ANCESTRY_INDEX_H
#define CURUPIXA_ANCESTRY_INDEX_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_ancestry_index crpx_ancestry_index;
/* Immutable snapshot of connectivity only: O(nodes log nodes) indices, O(log nodes)
 * LCA queries. No tree pointer, labels, lengths or all-pairs MRCA table is retained.
 * IDs refer to the input snapshot. Rebuild explicitly after changing connectivity. */
crpx_status crpx_ancestry_index_new(const crpx_tree *tree,crpx_ancestry_index **out);
void crpx_ancestry_index_free(crpx_ancestry_index *index);
size_t crpx_ancestry_index_nodes(const crpx_ancestry_index *index);
size_t crpx_ancestry_index_storage_size(const crpx_ancestry_index *index);
crpx_node_id crpx_ancestry_index_parent(const crpx_ancestry_index *index,crpx_node_id node);
uint32_t crpx_ancestry_index_depth(const crpx_ancestry_index *index,crpx_node_id node);
bool crpx_ancestry_index_is_leaf(const crpx_ancestry_index *index,crpx_node_id node);
bool crpx_ancestry_index_is_binary(const crpx_ancestry_index *index);
crpx_status crpx_ancestry_index_lca(const crpx_ancestry_index *index,crpx_node_id a,crpx_node_id b,
                                   crpx_node_id *result);
#ifdef __cplusplus
}
#endif
#endif
