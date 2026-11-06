/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose tree operations reimplemented from biomcmc-lib's topology and Newick modules.
 */
#ifndef CRPX_TREE_H
#define CRPX_TREE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef uint32_t crpx_node_id;
#define CRPX_NODE_NONE UINT32_MAX
typedef struct crpx_tree crpx_tree;
typedef crpx_tree *crpx_tree_t;
typedef struct crpx_tree_workspace crpx_tree_workspace;

/* An empty tree has no root. Node indices remain stable for the lifetime of a tree.
 * All allocation failures return a status; destruction accepts NULL.
 * Trees own their labels and branch lengths. No global context is required.
 * A read-only tree may be shared by threads; mutation and destruction require caller synchronisation.
 */
crpx_status crpx_tree_new(crpx_tree **result);
void crpx_tree_free(crpx_tree *tree);
crpx_status crpx_tree_clone(const crpx_tree *tree, crpx_tree **result);
crpx_status crpx_tree_reserve(crpx_tree *tree, size_t capacity);
/* The first node uses parent=CRPX_NODE_NONE. Later nodes are appended to a parent's child list.
 * Appending costs O(number of existing children); Newick construction remains O(number of nodes).
 */
crpx_status crpx_tree_add_node(crpx_tree *tree, crpx_node_id parent, crpx_node_id *result);
size_t crpx_tree_node_count(const crpx_tree *tree);
size_t crpx_tree_leaf_count(const crpx_tree *tree);
crpx_node_id crpx_tree_root(const crpx_tree *tree);
/* The relationship queries return CRPX_NODE_NONE for an absent relationship or invalid node. */
crpx_node_id crpx_tree_parent(const crpx_tree *tree, crpx_node_id node);
crpx_node_id crpx_tree_first_child(const crpx_tree *tree, crpx_node_id node);
crpx_node_id crpx_tree_next_sibling(const crpx_tree *tree, crpx_node_id node);
uint64_t crpx_tree_revision(const crpx_tree *tree);
/* Includes owned capacities and strings, excluding allocator overhead. */
size_t crpx_tree_storage_size(const crpx_tree *tree);
crpx_status crpx_tree_validate(const crpx_tree *tree);
/* Reparenting rejects cycles and moving the root. A node's branch length follows its incoming edge.
 * Rerooting at an internal node reverses path edges and shifts their lengths; the old root's stem length
 * follows the new root. A tip target or a unary old root is rejected to preserve the leaf/taxon set.
 * Rerooting may create a polytomy at the new root and a unary node at the previous root.
 */
crpx_status crpx_tree_reparent(crpx_tree *tree, crpx_node_id node, crpx_node_id parent);
crpx_status crpx_tree_reroot(crpx_tree *tree, crpx_node_id node);
/* Exchange disjoint non-root subtrees. Node-attached lengths/labels follow the exchanged nodes.
 * Siblings yield a no-op because unordered topology is unchanged.
 */
crpx_status crpx_tree_swap_subtrees(crpx_tree *tree, crpx_node_id first, crpx_node_id second);
/* Rooted binary SPR reuses prune's parent as the new graft node. regraft identifies the edge above
 * a node outside the pruned subtree; graft_fraction splits that edge between the graft and regraft.
 * Suppression adds the old parent/sister edge lengths, preserving paths amongst retained leaves.
 * Missing lengths propagate as missing. At a root graft, the new edge to the old root has length zero;
 * the stem follows the new root. Regrafting onto prune's parent/sister is rejected.
 * Internal labels remain node-attached; callers must reassess support values after rearrangement.
 */
crpx_status crpx_tree_spr(crpx_tree *tree, crpx_node_id prune, crpx_node_id regraft, double graft_fraction);
crpx_status crpx_tree_set_label(crpx_tree *tree, crpx_node_id node, const char *label);
const char *crpx_tree_label(const crpx_tree *tree, crpx_node_id node);
/* Branch lengths may be any finite value. Missing lengths return CRPX_ERROR_NOT_FOUND.
 * On failure, output arguments remain unchanged, except constructor outputs which become NULL.
 */
crpx_status crpx_tree_set_branch_length(crpx_tree *tree, crpx_node_id node, double length);
crpx_status crpx_tree_clear_branch_length(crpx_tree *tree, crpx_node_id node);
crpx_status crpx_tree_branch_length(const crpx_tree *tree, crpx_node_id node, double *length);

/* Parse exactly one semicolon-terminated Newick tree. Supports polytomies, unary nodes,
 * unnamed nodes, internal labels, nested bracket comments and single-quoted labels with doubled quotes.
 * Comments are skipped; underscores are preserved. Error offsets are zero-based byte positions.
 * Numbers use a decimal point irrespective of the process locale. Parsing/writing use iterative traversal.
 */
crpx_status crpx_tree_from_newick(const char *newick, crpx_tree **result, size_t *error_offset);
crpx_status crpx_tree_to_newick(const crpx_tree *tree, char **result);
void crpx_tree_string_free(char *string);

enum {
  CRPX_TREE_WORKSPACE_TRAVERSAL = 1u,
  CRPX_TREE_WORKSPACE_DEPTH = 2u,
  CRPX_TREE_WORKSPACE_SISTERS = 4u,
  CRPX_TREE_WORKSPACE_UPDATED = 8u,
  CRPX_TREE_WORKSPACE_ALL = 15u
};
enum { CRPX_TREE_UPDATED_UP = 1u, CRPX_TREE_UPDATED_DOWN = 2u, CRPX_TREE_UPDATED_BOTH = 3u };

/* Workspaces are separate, caller-owned analysis add-ons. Use one per concurrent analysis/thread.
 * Refresh computes selected arrays and clears update flags after any tree revision; a no-op refresh
 * preserves flags. Arrays are indexed by stable node ID, except preorder/postorder which contain IDs.
 * A sister is defined only for a node with exactly one sibling; polytomies yield CRPX_NODE_NONE.
 * Retain the bound tree until destroying or refreshing the workspace against a different live tree.
 */
crpx_status crpx_tree_workspace_new(uint32_t features, crpx_tree_workspace **result);
void crpx_tree_workspace_free(crpx_tree_workspace *workspace);
/* Release the borrowed tree identity while retaining allocated arrays for reuse. Call before freeing
 * a bound tree if the workspace must survive; it then reports zero nodes and requires refresh.
 */
void crpx_tree_workspace_reset(crpx_tree_workspace *workspace);
crpx_status crpx_tree_workspace_refresh(crpx_tree_workspace *workspace, const crpx_tree *tree);
bool crpx_tree_workspace_is_current(const crpx_tree_workspace *workspace, const crpx_tree *tree);
size_t crpx_tree_workspace_node_count(const crpx_tree_workspace *workspace);
const crpx_node_id *crpx_tree_workspace_preorder(const crpx_tree_workspace *workspace);
const crpx_node_id *crpx_tree_workspace_postorder(const crpx_tree_workspace *workspace);
const uint32_t *crpx_tree_workspace_depth(const crpx_tree_workspace *workspace);
const crpx_node_id *crpx_tree_workspace_sisters(const crpx_tree_workspace *workspace);
crpx_status crpx_tree_workspace_set_updated(crpx_tree_workspace *workspace, const crpx_tree *tree,
                                          crpx_node_id node, uint8_t flags);
crpx_status crpx_tree_workspace_get_updated(const crpx_tree_workspace *workspace, const crpx_tree *tree,
                                          crpx_node_id node, uint8_t *flags);

#ifdef __cplusplus
}
#endif
#endif
