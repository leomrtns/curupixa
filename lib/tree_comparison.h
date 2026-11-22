/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file tree_comparison.h
    \brief Prepare explicit tip correspondences for repeated split, length and rooted-vector comparisons.
    Robinson and Foulds (1981), doi:10.1016/0025-5564(81)90043-2;
    Kendall and Colijn (2016), doi:10.1093/molbev/msw124. */
#ifndef CURUPIXA_TREE_COMPARISON_H
#define CURUPIXA_TREE_COMPARISON_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct crpx_tree_comparison crpx_tree_comparison;
typedef struct {
  size_t taxa, first_tree_tips, second_tree_tips, first_splits, second_splits, rf;
  double weighted_rf, branch_score, length_scale, branch_scale;
} crpx_tree_distance_result;
/* tips lists distinct leaves in common taxon order. Other leaves remain in the
   input tree but do not contribute coordinates; this choice is explicit, never
   inferred by intersecting labels. Complete nonnegative lengths required when
   with_lengths=true. root stem ignored. Unary split lengths are summed; the two
   sides of a binary storage root are combined for unrooted comparisons. */
crpx_status crpx_tree_comparison_new (const crpx_tree *tree, const crpx_node_id *tips, size_t taxa, bool rooted,
                                      bool with_lengths, crpx_tree_comparison **out);
void crpx_tree_comparison_free (crpx_tree_comparison *comparison);
/* RF counts nontrivial unique splits; weighted RF sums absolute differences of
   edge lengths, including pendant edges. Branch score is the Euclidean norm of
   those differences, with absent splits assigned length zero. length_scale is
   the sum of edge lengths in both representations; branch_scale is the square root
   of the sum of squared lengths in both, for an optional branch-score denominator. No default normalisation:
   rf/(first_splits+second_splits) and weighted_rf/length_scale are available when
   denominators are positive. Length scores are NAN if lengths were not prepared. */
crpx_status crpx_tree_comparison_distance (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                                           crpx_tree_distance_result *result);
/* KC uses (1-lambda)*MRCA edge depth + lambda*MRCA root distance for each pair;
   pendant coordinates are (1-lambda)+lambda*pendant_length. Euclidean distance.
   It is a rooted metric, not a weighted sum of RF and path difference. */
crpx_status crpx_tree_comparison_kendall_colijn (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                                                 double lambda, double *distance);
/* Euclidean pairwise path difference; weighted=false counts edges, true uses lengths. */
crpx_status crpx_tree_comparison_path_difference (const crpx_tree_comparison *first, const crpx_tree_comparison *second,
                                                  bool weighted, double *distance);
/* One-off counterparts using the same complete-label policy as ordinary RF.
   Use prepared comparisons above for repeated scoring. No tip intersection. */
crpx_status crpx_tree_weighted_distances (const crpx_tree *first, const crpx_tree *second, bool rooted,
                                          crpx_tree_distance_result *result);
crpx_status crpx_tree_kendall_colijn_distance (const crpx_tree *first, const crpx_tree *second, double lambda,
                                               double *distance);
crpx_status crpx_tree_path_difference (const crpx_tree *first, const crpx_tree *second, bool weighted,
                                       double *distance);
#ifdef __cplusplus
}
#endif
#endif
