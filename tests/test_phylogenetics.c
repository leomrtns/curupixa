/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "parsimony.h"
#include "reconciliation.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK(expression) assert((expression) == CRPX_OK)

static crpx_node_id find_leaf(const crpx_tree *tree, const char *name)
{
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    const char *label = crpx_tree_label(tree, (crpx_node_id)i);
    if (label && !strcmp(label, name) && crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE)
      return (crpx_node_id)i;
  }
  return CRPX_NODE_NONE;
}

static void test_parsimony(void)
{
  crpx_tree *tree;
  crpx_parsimony_workspace *w;
  crpx_alignment *a;
  uint64_t score;
  OK(crpx_tree_from_newick("((A,B),(C,D));", &tree, NULL));
  OK(crpx_parsimony_workspace_new(&w));
  OK(crpx_alignment_from_text(">A\nAAAN\n>B\nAAGN\n>C\nCCAN\n>D\nCCGN\n", true, &a));
  OK(crpx_parsimony_score_alignment(tree, a, w, &score)); assert(score == 4);
  /* Reuse the analysis workspace after releasing its previous tree. */
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A,B,C,D);", &tree, NULL));
  OK(crpx_parsimony_score_alignment(tree, a, w, &score)); assert(score == 6);
  crpx_node_id tips[] = {find_leaf(tree,"A"), find_leaf(tree,"B"), find_leaf(tree,"C"), find_leaf(tree,"D")};
  uint64_t masks[] = {1, 2, 3, 3}, per;
  OK(crpx_parsimony_score(tree, tips, 4, 2, 1, masks, NULL, w, &score, &per)); assert(score == 1 && per == 1);
  masks[0] = 4;
  score = 77; per = 88;
  assert(crpx_parsimony_score(tree, tips, 4, 2, 1, masks, NULL, w, &score, &per) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(score == 77 && per == 88);
#if SIZE_MAX == UINT64_MAX
  uint64_t overflow_masks[] = {1, 2, 3, 3, 1, 2, 3, 3}, pattern_scores[2] = {88, 99};
  size_t overflow_weights[] = {SIZE_MAX, 1};
  assert(crpx_parsimony_score(tree, tips, 4, 2, 2, overflow_masks, overflow_weights, w, &score, pattern_scores)
         == CRPX_ERROR_OVERFLOW);
  assert(score == 77 && pattern_scores[0] == 88 && pattern_scores[1] == 99);
#endif
  crpx_parsimony_workspace_free(w); crpx_tree_free(tree); crpx_alignment_free(a);
}

static void reconcile_case(const char *newick, uint64_t duplications, uint64_t losses, uint64_t extra)
{
  crpx_tree *gene, *species;
  crpx_reconciliation *r;
  OK(crpx_tree_from_newick(newick, &gene, NULL));
  OK(crpx_tree_from_newick("((A,B),C);", &species, NULL));
  size_t nodes = crpx_tree_node_count(gene);
  crpx_node_id *mapping = malloc(nodes*sizeof(*mapping)); assert(mapping);
  for (size_t i = 0; i < nodes; ++i) {
    const char *label = crpx_tree_label(gene, (crpx_node_id)i);
    char species_name[2] = {label && *label ? *label : '?', 0};
    mapping[i] = find_leaf(species, species_name);
  }
  OK(crpx_reconcile_lca(gene, species, mapping, &r));
  assert(crpx_reconciliation_duplications(r) == duplications);
  assert(crpx_reconciliation_losses(r) == losses);
  assert(crpx_reconciliation_extra_lineages(r) == extra);
  crpx_reconciliation_free(r); free(mapping); crpx_tree_free(gene); crpx_tree_free(species);
}

int main(void)
{
  test_parsimony();
  reconcile_case("((A,B),C);", 0, 0, 0);
  reconcile_case("((A1,A2),B);", 1, 0, 0);
  reconcile_case("((A,C),B);", 1, 3, 1);
  reconcile_case("(A,C);", 0, 1, 0);
  puts("parsimony and rooted reconciliation: passed");
  return 0;
}
