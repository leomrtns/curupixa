/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * Reimplementation of biomcmc-lib Newick/tree collections with explicit ownership.
 */
#ifndef CRPX_NEWICK_SPACE_H
#define CRPX_NEWICK_SPACE_H
#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_tree_space crpx_tree_space;
/* Collections own independent tree copies. Concurrent readers may share a const collection;
 * mutation requires caller synchronisation. Returned views live until collection mutation/free.
 */
crpx_status crpx_tree_space_new(crpx_tree_space **result);
void crpx_tree_space_free(crpx_tree_space *space);
crpx_status crpx_tree_space_clone(const crpx_tree_space *space, crpx_tree_space **result);
size_t crpx_tree_space_size(const crpx_tree_space *space);
uint64_t crpx_tree_space_total_weight(const crpx_tree_space *space);
const crpx_tree *crpx_tree_space_get(const crpx_tree_space *space, size_t index);
crpx_status crpx_tree_space_weight(const crpx_tree_space *space, size_t index, uint64_t *weight);
crpx_status crpx_tree_space_add(crpx_tree_space *space, const crpx_tree *tree, uint64_t weight);
/* Appending a list is transactional: errors preserve the original collection.
 * Empty input/trailing whitespace/comments are allowed; each tree requires a semicolon.
 */
crpx_status crpx_tree_space_add_newick(crpx_tree_space *space, const char *text, size_t *error_offset);
/* Auto-detect Newick or NEXUS (#NEXUS, TREES blocks, TREE/UTREE and TRANSLATE).
 * NEXUS tree names/root annotations are not retained; the backbone stores the supplied rooting.
 */
crpx_status crpx_tree_space_add_trees(crpx_tree_space *space, const char *text, size_t *error_offset);
/* File reading detects Newick/NEXUS and compression through crpx_file_open(AUTO), streaming one tree at a time.
 * A failed read destroys the partial result. Offsets count decompressed bytes from the start.
 */
crpx_status crpx_tree_space_read_file(const char *path, crpx_tree_space **result, size_t *error_offset);
/* Writes one representative per entry; weights are available separately and are not encoded in Newick.
 * Compression follows the filename suffix, subject to available optional compression libraries.
 */
crpx_status crpx_tree_space_write_file(const crpx_tree_space *space, const char *path);
/* Merge topologically identical entries, retaining the first tree's labels/lengths and summing weights.
 * Requires unique, matching leaf-label sets. Rooted/unrooted comparison uses exact RF distance.
 * Transactional and O(entries^2) RF comparisons; use for modest collections, not an index for huge samples.
 */
crpx_status crpx_tree_space_compact(crpx_tree_space *space, bool rooted);
/* Select original entry indices burnin, burnin+thin, ... before compaction. thin must be positive.
 * Entry weights are retained; selection counts entries rather than expanding weighted repetitions.
 */
crpx_status crpx_tree_space_select(const crpx_tree_space *space, size_t burnin, size_t thin,
                                 crpx_tree_space **result);
/* Transactionally append independent copies, multiplying integer weights. Self-merge is supported. */
crpx_status crpx_tree_space_merge(crpx_tree_space *space, const crpx_tree_space *from, uint64_t multiplier);
crpx_status crpx_tree_space_frequency(const crpx_tree_space *space, size_t index, double *frequency);

#ifdef __cplusplus
}
#endif
#endif
