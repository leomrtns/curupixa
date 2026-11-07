/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose replacements for biomcmc-lib topology randomisation.
 */
#ifndef CRPX_TOPOLOGY_RANDOMISE_H
#define CRPX_TOPOLOGY_RANDOMISE_H
#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The caller owns the generator/state, so MPI ranks and threads select their own streams.
 * The callback must supply uniform 64-bit words; bounded draws use rejection sampling.
 */
typedef uint64_t (*crpx_tree_random_function)(void *state);
/* Randomly coalesce pairs of active lineages. Uniform pair choices give unequal probabilities
 * to labelled rooted trees. Leaves receive supplied labels or
 * zero-based decimal names. Branch lengths remain absent. Construction costs O(leaves).
 */
crpx_status crpx_tree_random_binary(size_t leaves, const char *const *labels,
                                   crpx_tree_random_function random_word, void *state, crpx_tree **result);
/* Rooted binary NNI exchanges the selected child (0=first, 1=second) with its parent's other child.
 * The selected node and its parent must both have exactly two children.
 */
crpx_status crpx_tree_nni(crpx_tree *tree, crpx_node_id node, unsigned child_choice);
/* Uniformly choose an eligible internal edge, then one of its two rooted NNI alternatives. */
crpx_status crpx_tree_random_nni(crpx_tree *tree, crpx_tree_random_function random_word, void *state);
/* Uniformly choose an eligible prune node, then a legal regraft edge for that prune.
 * This proposal is not uniform over all distinct SPR neighbours. The branch-length policy is
 * crpx_tree_spr's explicit graft_fraction policy. Selection costs O(nodes) time and memory.
 */
crpx_status crpx_tree_random_spr(crpx_tree *tree, double graft_fraction,
                               crpx_tree_random_function random_word, void *state);

#ifdef __cplusplus
}
#endif
#endif
