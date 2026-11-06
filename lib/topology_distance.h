/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose replacements for biomcmc-lib tree distance operations.
 */
#ifndef CRPX_TOPOLOGY_DISTANCE_H
#define CRPX_TOPOLOGY_DISTANCE_H
#include "tree.h"
#include "distance_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Ancestry includes the node itself. Queries cost O(tree height), without allocating a workspace. */
crpx_status crpx_tree_is_ancestor(const crpx_tree *tree, crpx_node_id ancestor, crpx_node_id node, bool *result);
crpx_status crpx_tree_lca(const crpx_tree *tree, crpx_node_id first, crpx_node_id second, crpx_node_id *result);
/* The root's stem length is excluded. NaN missing_length requests an error for absent edges;
 * a finite missing_length supplies a default. Output arguments remain unchanged on failure.
 */
crpx_status crpx_tree_patristic_distance(const crpx_tree *tree, crpx_node_id first, crpx_node_id second,
                                       double missing_length, double *result);
/* Leaves use increasing node-ID order. The output is packed symmetric storage.
 * This calculation requires non-negative edge lengths/defaults and costs O(leaves * nodes).
 * The optional leaf_order output must be released with crpx_tree_node_ids_free().
 */
crpx_status crpx_tree_patristic_matrix(const crpx_tree *tree, double missing_length,
                                     crpx_distance_matrix_t *result, crpx_node_id **leaf_order);
void crpx_tree_node_ids_free(crpx_node_id *nodes);
/* Exact Robinson-Foulds distance: count clades/splits present in exactly one tree.
 * Leaves must have unique labels, and both trees must contain the same labels.
 * Internal labels, child order, branch lengths, unary duplicate splits and trivial splits are ignored.
 * rooted=true compares rooted clades; false compares unrooted bipartitions.
 * Separate temporary bitsets cost O(nodes * ceil(leaves/64)) storage; the backbone is unchanged.
 */
crpx_status crpx_tree_rf_distance(const crpx_tree *first, const crpx_tree *second, bool rooted, size_t *result);

#ifdef __cplusplus
}
#endif
#endif
