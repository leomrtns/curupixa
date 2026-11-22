/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file rooting.h
    \brief Locate midpoint, minimum-variance and minimal ancestor deviation roots without changing a tree.
 */
#ifndef CURUPIXA_ROOTING_H
#define CURUPIXA_ROOTING_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_rooting_workspace crpx_rooting_workspace;
typedef struct {
  crpx_node_id child;
  double fraction, score;
} crpx_root_location;
/* Reusable O(nodes) scratch; independent workspaces permit concurrent calls on a read-only tree.
   No tree identity is retained. NULL workspace in a calculation uses temporary storage. */
crpx_status crpx_rooting_workspace_new (crpx_rooting_workspace **out);
void crpx_rooting_workspace_free (crpx_rooting_workspace *workspace);
/* All methods require a valid tree with at least two tips and complete finite nonnegative edge lengths.
   The root stem is ignored. Polytomies and non-root unary vertices are supported; a unary storage root is rejected
   because it would introduce a non-taxon terminal vertex. Existing binary storage roots merely subdivide an edge.
   child identifies an input edge, fraction is measured from its parent in [0,1]. Results are invalidated by tree
   edits. Apply with crpx_tree_reroot_edge(), choosing the old-root policy explicitly. Input is never modified and
   output is unchanged on failure. Equal scores choose the lowest child ID (then the smaller fraction); tied vertices
   can be represented on several edges. No biological rooting flag or topology constraint is changed.
   Lengths need not be ultrametric. OVERFLOW includes an unrepresentable score. */
/* Midpoint of a diameter path. score is the maximum root-to-tip distance, in input length units. O(nodes) time.
   Mai, Sayyari & Mirarab (2017), doi:10.1371/journal.pone.0182238. */
crpx_status crpx_tree_root_midpoint (const crpx_tree *tree, crpx_rooting_workspace *workspace,
                                     crpx_root_location *out);
/* Globally minimise population variance of root-to-tip distances. score has squared input length units.
   O(nodes) time. Mai, Sayyari & Mirarab (2017), doi:10.1371/journal.pone.0182238. */
crpx_status crpx_tree_root_minvar (const crpx_tree *tree, crpx_rooting_workspace *workspace,
                                   crpx_root_location *out);
/* Globally minimise RMS relative pairwise ancestor deviation (MAD A1-B1). O(nodes*tips) time and O(nodes) scratch.
   Each pair is scored at its LCA induced by the candidate root, including pairs on the same side of its edge.
   Zero-distance tip pairs are omitted from both sum and denominator; an entirely zero-distance tree is invalid.
   This explicit zero-path convention does not perturb edge lengths. score is dimensionless, in [0,1].
   Tria, Landan & Dagan (2017), doi:10.1038/s41559-017-0193. */
crpx_status crpx_tree_root_mad (const crpx_tree *tree, crpx_rooting_workspace *workspace,
                                crpx_root_location *out);
#ifdef __cplusplus
}
#endif
#endif
