/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file sampling.h
    \brief Reusable draws and tree construction, independent of biological prior models. */
#ifndef CURUPIXA_SAMPLING_H
#define CURUPIXA_SAMPLING_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef double (*crpx_uniform_function) (void *state);
/* Callback returns a finite uniform in [0,1). Each sampling operation checks it. */
crpx_status crpx_sample_uniform (crpx_uniform_function uniform, void *state, double *value);
crpx_status crpx_sample_log_weights (const double *weights, size_t count, crpx_uniform_function uniform, void *state,
                                     size_t *selected);
crpx_status crpx_sample_permutation (size_t *indices, size_t count, crpx_uniform_function uniform, void *state);
/* Import a parent array in O(nodes) time. Exactly one parent must be NODE_NONE;
   remaining entries identify input node indices. Child order follows input index.
   Labels and lengths are node-indexed; NULL omits them, NAN marks missing lengths.
   The root's optional length is its stem. Input indices need not be topologically
   sorted; output IDs may differ. Internal labels and unary nodes are supported. */
crpx_status crpx_tree_from_parents (const crpx_node_id *parents, size_t nodes, const char *const *labels,
                                    const double *lengths, crpx_tree **out);
/* Leaves are [0,leaves), internal vertices [leaves,2*leaves-1). children has two
   entries per internal vertex. Every non-root vertex must occur exactly once.
   ages=NULL omits lengths, otherwise ages are finite and parent>=child. labels
   is leaf-indexed and may be NULL. Node IDs in the resulting tree are independent
   of these input indices. No input pointers survive construction. */
crpx_status crpx_tree_from_binary_history (size_t leaves, const size_t *children, size_t root, const double *ages,
                                           const char *const *labels, crpx_tree **out);
/* Uniformly resolve each multifurcation over labelled rooted binary resolutions.
   Preserve each existing incoming length/label; inserted edges have length zero.
   Unary nodes are retained. Uniformity is conditional on the supplied clades. */
crpx_status crpx_tree_resolve_random (const crpx_tree *tree, crpx_uniform_function uniform, void *state,
                                      crpx_tree **out);
/* Insert a new labelled tip on edge child's incoming edge; child=root inserts
   above the root. Non-root edge length is split by fraction measured from parent.
   At a root insertion stem_length supplies the new edge to the old root. The
   pendant length is supplied separately. Input is unchanged; output node IDs may
   differ. All existing non-root lengths must be finite and nonnegative. */
crpx_status crpx_tree_insert_tip (const crpx_tree *tree, crpx_node_id child, double fraction, double pendant_length,
                                  double stem_length, const char *label, crpx_tree **out);
#ifdef __cplusplus
}
#endif
#endif
