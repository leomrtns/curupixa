/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "lowlevel_algorithms.h"
#include "splitset_distances.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

static double brute_assignment(const double *costs, size_t rows, size_t columns, size_t row, unsigned int used)
{
  if (row == rows) return 0;
  double best = INFINITY;
  for (size_t j = 0; j < columns; ++j) if (!(used & (1u << j))) {
    double cost = costs[row*columns+j]+brute_assignment(costs, rows, columns, row+1, used | (1u << j));
    if (cost < best) best = cost;
  }
  return best;
}
static void test_assignment(void)
{
  const double largest = DBL_MAX;
  double extreme_total;
  size_t extreme_assignment = 99;
  assert(crpx_assignment_solve(&largest, 1, 1, &extreme_assignment, &extreme_total) == CRPX_OK);
  assert(extreme_total == DBL_MAX && extreme_assignment == 0);
  for (size_t rows = 1; rows <= 5; ++rows) for (size_t cols = rows; cols <= 6; ++cols) {
    double costs[30], total;
    size_t assignment[5];
    for (size_t i = 0; i < rows*cols; ++i) costs[i] = (double)((i*37+rows*13+i*i*19) % 23);
    assert(crpx_assignment_solve(costs, rows, cols, assignment, &total) == CRPX_OK);
    assert(total == brute_assignment(costs, rows, cols, 0, 0));
    unsigned int used = 0;
    for (size_t i = 0; i < rows; ++i) { assert(!(used & (1u << assignment[i]))); used |= 1u << assignment[i]; }
  }
}
static void test_edit_distance(void)
{
  uint64_t distance;
  assert(crpx_levenshtein_distance("kitten", 6, "sitting", 7, 1, 1, false, &distance) == CRPX_OK && distance == 3);
  assert(crpx_levenshtein_distance("ACGTT", 5, "ACG", 3, 1, 1, true, &distance) == CRPX_OK && distance == 0);
  assert(crpx_levenshtein_distance("A", 1, "G", 1, 1, 1, true, &distance) == CRPX_OK && distance == 1);
  assert(crpx_levenshtein_distance("A", 1, "G", 1, UINT64_MAX, UINT64_MAX, false, &distance) == CRPX_OK);
  assert(distance == UINT64_MAX);
  assert(crpx_levenshtein_distance("AA", 2, "GG", 2, UINT64_MAX, UINT64_MAX, false, &distance) == CRPX_ERROR_OVERFLOW);
  assert(crpx_levenshtein_distance("A", 1, "A", 1, UINT64_MAX, UINT64_MAX, false, &distance) == CRPX_OK && !distance);
  assert(crpx_levenshtein_distance(NULL, 0, NULL, 0, 1, 1, false, &distance) == CRPX_OK && !distance);
  const char *repeated = "aaaa";
  assert(crpx_string_last_substring(repeated, "aa") == repeated+2);
  assert(crpx_string_common_prefix("ACGT", "ACGG") == 3);
}
static crpx_node_id leaf(const crpx_tree *tree, const char *label)
{
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    const char *s = crpx_tree_label(tree, (crpx_node_id)i);
    if (s && !strcmp(s, label)) return (crpx_node_id)i;
  }
  return CRPX_NODE_NONE;
}
static void test_splits(void)
{
  crpx_tree *a, *b, *star;
  crpx_splitset *x, *y, *z;
  assert(crpx_tree_from_newick("((A,B),(C,D));", &a, NULL) == CRPX_OK);
  assert(crpx_tree_from_newick("((A,C),(B,D));", &b, NULL) == CRPX_OK);
  assert(crpx_tree_from_newick("(A,B,C,D);", &star, NULL) == CRPX_OK);
  const char *names[] = {"A", "B", "C", "D"};
  crpx_node_id tips[4];
  for (size_t i = 0; i < 4; ++i) tips[i] = leaf(a, names[i]);
  assert(crpx_splitset_from_tree(a, tips, 4, &x) == CRPX_OK);
  for (size_t i = 0; i < 4; ++i) tips[i] = leaf(b, names[i]);
  assert(crpx_splitset_from_tree(b, tips, 4, &y) == CRPX_OK);
  for (size_t i = 0; i < 4; ++i) tips[i] = leaf(star, names[i]);
  assert(crpx_splitset_from_tree(star, tips, 4, &z) == CRPX_OK);
  size_t rf;
  uint64_t h;
  assert(crpx_splitset_compare(x, y, &rf, &h) == CRPX_OK && rf == 2 && h == 2);
  assert(crpx_splitset_compare(x, z, &rf, &h) == CRPX_OK && rf == 1 && h == 2);
  assert(crpx_splitset_compare(x, x, &rf, &h) == CRPX_OK && rf == 0 && h == 0);
  crpx_splitset_free(x); crpx_splitset_free(y); crpx_splitset_free(z);
  crpx_tree_free(b); crpx_tree_free(star);
  assert(crpx_tree_from_newick("(((A1,A2),B),(C,D));", &b, NULL) == CRPX_OK);
  crpx_node_id mapping[9];
  assert(crpx_tree_node_count(b) == 9);
  for (size_t node = 0; node < 9; ++node) {
    const char *label = crpx_tree_label(b, (crpx_node_id)node);
    char name[2] = {label && *label ? *label : '?', 0}; mapping[node] = leaf(a, name);
  }
  assert(crpx_gene_species_split_distances(b, a, mapping, &rf, &h) == CRPX_OK && rf == 0 && h == 0);
  crpx_tree_free(a); crpx_tree_free(b);
}
static void balanced_subtree(crpx_tree *tree, crpx_node_id parent, crpx_node_id *tips, size_t first, size_t count)
{
  crpx_node_id node;
  assert(crpx_tree_add_node(tree, parent, &node) == CRPX_OK);
  if (count == 1) { tips[first] = node; return; }
  balanced_subtree(tree, node, tips, first, count/2);
  balanced_subtree(tree, node, tips, first+count/2, count-count/2);
}
static void test_large_rf(void)
{
  enum { TAXA = 128 };
  crpx_tree *tree, *other;
  crpx_node_id tips[TAXA];
  crpx_splitset *a, *b;
  assert(crpx_tree_new(&tree) == CRPX_OK);
  balanced_subtree(tree, CRPX_NODE_NONE, tips, 0, TAXA);
  assert(crpx_tree_clone(tree, &other) == CRPX_OK);
  assert(crpx_splitset_from_tree(tree, tips, TAXA, &a) == CRPX_OK);
  /* The two root edges are complementary; extraction must deduplicate them. */
  assert(crpx_splitset_count(a) == TAXA-3);
  assert(crpx_tree_swap_subtrees(other, tips[0], tips[TAXA-1]) == CRPX_OK);
  assert(crpx_splitset_from_tree(other, tips, TAXA, &b) == CRPX_OK);
  size_t identical = SIZE_MAX, forward = 0, reverse = 0;
  assert(crpx_splitset_compare(a, a, &identical, NULL) == CRPX_OK && !identical);
  assert(crpx_splitset_compare(a, b, &forward, NULL) == CRPX_OK && forward);
  assert(crpx_splitset_compare(b, a, &reverse, NULL) == CRPX_OK && reverse == forward);
  crpx_splitset_free(a); crpx_splitset_free(b); crpx_tree_free(tree); crpx_tree_free(other);
}
int main(void)
{
  test_assignment(); test_edit_distance(); test_splits(); test_large_rf();
  puts("assignment, edit and mapped split distances: passed");
  return 0;
}
