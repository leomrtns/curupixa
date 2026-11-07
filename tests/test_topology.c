/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 */
#include "topology_distance.h"
#include "topology_randomise.h"
#include "upgma.h"
#include "newick_space.h"
#include "file_compression.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(test) do { if (!(test)) { \
  fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #test); abort(); \
} } while (0)
#define OK(expression) CHECK((expression) == CRPX_OK)

static uint64_t
random_word(void *state)
{
  uint64_t *word = (uint64_t *) state;
  *word += UINT64_C(0x9e3779b97f4a7c15);
  {
    uint64_t mixed = *word;
    mixed = (mixed ^ (mixed >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    mixed = (mixed ^ (mixed >> 27)) * UINT64_C(0x94d049bb133111eb);
    return mixed ^ (mixed >> 31);
  }
}

static crpx_node_id
find_label(const crpx_tree *tree, const char *label)
{
  size_t i;
  for (i = 0; i < crpx_tree_node_count(tree); ++i) {
    const char *name = crpx_tree_label(tree, (crpx_node_id) i);
    if (name && !strcmp(name, label)) return (crpx_node_id) i;
  }
  return CRPX_NODE_NONE;
}

static void
check_rf(const char *first, const char *second, size_t rooted, size_t unrooted)
{
  crpx_tree *a, *b;
  size_t distance;
  OK(crpx_tree_from_newick(first, &a, NULL));
  OK(crpx_tree_from_newick(second, &b, NULL));
  OK(crpx_tree_rf_distance(a, b, true, &distance)); CHECK(distance == rooted);
  OK(crpx_tree_rf_distance(b, a, true, &distance)); CHECK(distance == rooted);
  OK(crpx_tree_rf_distance(a, b, false, &distance)); CHECK(distance == unrooted);
  OK(crpx_tree_rf_distance(b, a, false, &distance)); CHECK(distance == unrooted);
  crpx_tree_free(a); crpx_tree_free(b);
}

static void
test_rf(void)
{
  crpx_tree *a, *b, *other;
  size_t distance = 19, i;
  uint64_t seed = 123;
  check_rf("((A,B),(C,D));", "((D,C),(B,A));", 0, 0);
  check_rf("((A,B),(C,D));", "((A,C),(B,D));", 4, 2);
  check_rf("((A,B),(C,D));", "(A,B,C,D);", 2, 1);
  check_rf("((A,B),C);", "(((A,B)),C);", 0, 0);
  check_rf("((A,B),C);", "(A,(B,C));", 2, 0);
  check_rf("A;", "A;", 0, 0);
  OK(crpx_tree_from_newick("(A,A);", &a, NULL));
  OK(crpx_tree_from_newick("(A,B);", &b, NULL));
  CHECK(crpx_tree_rf_distance(a, b, false, &distance) == CRPX_ERROR_INVALID_ARGUMENT && distance == 19);
  crpx_tree_free(a);
  OK(crpx_tree_from_newick("(A,C);", &a, NULL));
  CHECK(crpx_tree_rf_distance(a, b, true, &distance) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(a); crpx_tree_free(b);
  OK(crpx_tree_random_binary(71, NULL, random_word, &seed, &a));
  OK(crpx_tree_clone(a, &b));
  OK(crpx_tree_rf_distance(a, b, false, &distance)); CHECK(!distance);
  for (i = 0; i < crpx_tree_node_count(b); ++i)
    if (i != crpx_tree_root(b) && crpx_tree_first_child(b, (crpx_node_id) i) != CRPX_NODE_NONE) break;
  OK(crpx_tree_reroot(b, (crpx_node_id) i));
  OK(crpx_tree_rf_distance(a, b, false, &distance)); CHECK(!distance);
  OK(crpx_tree_random_binary(71, NULL, random_word, &seed, &other));
  OK(crpx_tree_rf_distance(a, other, false, &distance)); CHECK(distance <= 2 * (71 - 3));
  OK(crpx_tree_rf_distance(a, other, true, &distance)); CHECK(distance <= 2 * (71 - 2));
  crpx_tree_free(a); crpx_tree_free(b); crpx_tree_free(other);
}

static void
test_distances(void)
{
  crpx_tree *tree;
  crpx_distance_matrix_t matrix;
  crpx_node_id a, b, c, x, root, lca, *leaves;
  size_t i, j;
  double distance, value;
  bool ancestor;
  OK(crpx_tree_from_newick("((A:1,B:2)X:3,C:4)R:100;", &tree, NULL));
  a = find_label(tree, "A"); b = find_label(tree, "B");
  c = find_label(tree, "C"); x = find_label(tree, "X"); root = crpx_tree_root(tree);
  OK(crpx_tree_lca(tree, a, b, &lca)); CHECK(lca == x);
  OK(crpx_tree_lca(tree, b, c, &lca)); CHECK(lca == root);
  OK(crpx_tree_lca(tree, x, a, &lca)); CHECK(lca == x);
  OK(crpx_tree_lca(tree, a, a, &lca)); CHECK(lca == a);
  OK(crpx_tree_is_ancestor(tree, x, a, &ancestor)); CHECK(ancestor);
  OK(crpx_tree_is_ancestor(tree, x, c, &ancestor)); CHECK(!ancestor);
  OK(crpx_tree_is_ancestor(tree, a, a, &ancestor)); CHECK(ancestor);
  OK(crpx_tree_patristic_distance(tree, a, b, NAN, &distance)); CHECK(distance == 3);
  OK(crpx_tree_patristic_distance(tree, a, c, NAN, &distance)); CHECK(distance == 8);
  OK(crpx_tree_patristic_distance(tree, root, c, NAN, &distance)); CHECK(distance == 4);
  OK(crpx_tree_patristic_matrix(tree, NAN, &matrix, &leaves));
  CHECK(crpx_distance_matrix_size(matrix) == 3);
  for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j) {
    OK(crpx_distance_matrix_get(matrix, i, j, &value));
    OK(crpx_tree_patristic_distance(tree, leaves[i], leaves[j], NAN, &distance));
    CHECK(distance == value);
  }
  crpx_distance_matrix_free(matrix); crpx_tree_node_ids_free(leaves);
  OK(crpx_tree_clear_branch_length(tree, a));
  CHECK(crpx_tree_patristic_distance(tree, a, b, NAN, &distance) == CRPX_ERROR_NOT_FOUND);
  OK(crpx_tree_patristic_distance(tree, a, b, 1.5, &distance)); CHECK(distance == 3.5);
  CHECK(crpx_tree_patristic_matrix(tree, NAN, &matrix, &leaves) == CRPX_ERROR_NOT_FOUND);
  CHECK(!matrix && !leaves);
  OK(crpx_tree_patristic_matrix(tree, 1, &matrix, NULL));
  crpx_distance_matrix_free(matrix);
  OK(crpx_tree_set_branch_length(tree, a, -1));
  CHECK(crpx_tree_patristic_matrix(tree, 1, &matrix, NULL) == CRPX_ERROR_INVALID_ARGUMENT && !matrix);
  crpx_tree_free(tree);
}

static void
test_spr_exhaustive(void)
{
  crpx_tree *original, *tree;
  crpx_node_id prune, target, a, b;
  size_t count, valid = 0, invalid = 0;
  bool below_a, below_b;
  crpx_status status;
  double before, after;
  char *old_text, *text;
  OK(crpx_tree_from_newick("(((A:1,B:2)X:3,C:4)Y:5,(D:6,E:7)Z:8)R:9;", &original, NULL));
  count = crpx_tree_node_count(original);
  OK(crpx_tree_to_newick(original, &old_text));
  for (prune = 0; prune < count; ++prune) for (target = 0; target < count; ++target) {
    OK(crpx_tree_clone(original, &tree));
    status = crpx_tree_spr(tree, prune, target, 0.25);
    OK(crpx_tree_validate(tree));
    if (status == CRPX_OK) {
      ++valid;
      CHECK(crpx_tree_node_count(tree) == count);
      CHECK(crpx_tree_leaf_count(tree) == 5);
      for (a = 0; a < count; ++a) {
        if (crpx_tree_first_child(original, a) != CRPX_NODE_NONE) continue;
        OK(crpx_tree_is_ancestor(original, prune, a, &below_a));
        if (below_a) continue;
        for (b = a + 1; b < count; ++b) {
          if (crpx_tree_first_child(original, b) != CRPX_NODE_NONE) continue;
          OK(crpx_tree_is_ancestor(original, prune, b, &below_b));
          if (below_b) continue;
          OK(crpx_tree_patristic_distance(original, a, b, NAN, &before));
          OK(crpx_tree_patristic_distance(tree, a, b, NAN, &after));
          CHECK(fabs(before - after) < 1e-10);
        }
      }
    }
    else {
      ++invalid;
      CHECK(status == CRPX_ERROR_TOPOLOGY);
      OK(crpx_tree_to_newick(tree, &text));
      CHECK(!strcmp(old_text, text));
      crpx_tree_string_free(text);
    }
    crpx_tree_free(tree);
  }
  CHECK(valid && invalid);
  crpx_tree_string_free(old_text);
  crpx_tree_free(original);
}

static void
test_randomisation(void)
{
  uint64_t seed = 11, same_seed = 11;
  crpx_tree *tree, *copy, *single;
  size_t i, distance;
  crpx_node_id node;
  char *first, *second;
  const char *labels[] = { "alpha", "beta", "gamma" };
  OK(crpx_tree_random_binary(64, NULL, random_word, &seed, &tree));
  OK(crpx_tree_random_binary(64, NULL, random_word, &same_seed, &copy));
  OK(crpx_tree_to_newick(tree, &first));
  OK(crpx_tree_to_newick(copy, &second));
  CHECK(!strcmp(first, second));
  crpx_tree_string_free(first); crpx_tree_string_free(second);
  CHECK(crpx_tree_node_count(tree) == 127 && crpx_tree_leaf_count(tree) == 64);
  for (i = 0; i < 127; ++i) if (i != crpx_tree_root(tree) && crpx_tree_first_child(tree, (crpx_node_id) i) != CRPX_NODE_NONE) break;
  node = (crpx_node_id) i;
  OK(crpx_tree_nni(tree, node, 0));
  OK(crpx_tree_validate(tree));
  OK(crpx_tree_rf_distance(tree, copy, true, &distance)); CHECK(distance == 2);
  CHECK(crpx_tree_nni(tree, crpx_tree_root(tree), 0) == CRPX_ERROR_TOPOLOGY);
  for (i = 0; i < 200; ++i) {
    OK(crpx_tree_random_nni(tree, random_word, &seed));
    OK(crpx_tree_random_spr(tree, 0.5, random_word, &seed));
    OK(crpx_tree_validate(tree));
    CHECK(crpx_tree_node_count(tree) == 127 && crpx_tree_leaf_count(tree) == 64);
  }
  OK(crpx_tree_rf_distance(tree, copy, false, &distance)); CHECK(distance <= 2 * (64 - 3));
  crpx_tree_free(tree); crpx_tree_free(copy);
  OK(crpx_tree_random_binary(3, labels, random_word, &seed, &tree));
  CHECK(find_label(tree, "alpha") != CRPX_NODE_NONE && find_label(tree, "gamma") != CRPX_NODE_NONE);
  crpx_tree_free(tree);
  OK(crpx_tree_random_binary(1, NULL, random_word, &seed, &single));
  CHECK(crpx_tree_random_nni(single, random_word, &seed) == CRPX_ERROR_NOT_FOUND);
  CHECK(crpx_tree_random_spr(single, 0.5, random_word, &seed) == CRPX_ERROR_NOT_FOUND);
  crpx_tree_free(single);
  CHECK(crpx_tree_random_binary(0, NULL, random_word, &seed, &tree) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
}

static void
compare_tree_matrix(const crpx_tree *tree, const struct crpx_distance_matrix_struct *matrix,
                    const char *const *labels, size_t count)
{
  size_t i, j;
  crpx_node_id a, b;
  double expected, observed;
  for (i = 0; i < count; ++i) for (j = 0; j < i; ++j) {
    a = find_label(tree, labels[i]); b = find_label(tree, labels[j]);
    CHECK(a != CRPX_NODE_NONE && b != CRPX_NODE_NONE);
    OK(crpx_distance_matrix_get(matrix, i, j, &expected));
    OK(crpx_tree_patristic_distance(tree, a, b, NAN, &observed));
    CHECK(fabs(expected - observed) < 1e-10);
  }
}

static void
test_distance_construction(void)
{
  crpx_tree *source, *tree, *other;
  crpx_distance_matrix_t distances;
  crpx_node_id *order;
  const char *labels[5];
  size_t i, distance;
  double length;
  bool negative;
  OK(crpx_tree_from_newick("((A:1,B:1):2,(C:2,D:2):1);", &source, NULL));
  OK(crpx_tree_patristic_matrix(source, NAN, &distances, &order));
  for (i = 0; i < 4; ++i) labels[i] = crpx_tree_label(source, order[i]);
  OK(crpx_tree_upgma(distances, labels, &tree));
  OK(crpx_tree_validate(tree));
  compare_tree_matrix(tree, distances, labels, 4);
  OK(crpx_tree_rf_distance(tree, source, true, &distance)); CHECK(!distance);
  crpx_tree_free(tree);
  OK(crpx_tree_single_linkage(distances, labels, &tree));
  compare_tree_matrix(tree, distances, labels, 4);
  crpx_tree_free(tree); crpx_tree_node_ids_free(order);
  crpx_tree_free(source); crpx_distance_matrix_free(distances);
  OK(crpx_tree_from_newick("((A:1,B:2):3,(C:4,(D:5,E:6):7):8);", &source, NULL));
  OK(crpx_tree_patristic_matrix(source, NAN, &distances, &order));
  for (i = 0; i < 5; ++i) labels[i] = crpx_tree_label(source, order[i]);
  OK(crpx_tree_neighbour_joining(distances, labels, true, &tree));
  OK(crpx_tree_bionj(distances, labels, true, &other));
  compare_tree_matrix(tree, distances, labels, 5);
  compare_tree_matrix(other, distances, labels, 5);
  OK(crpx_tree_rf_distance(tree, source, false, &distance)); CHECK(!distance);
  OK(crpx_tree_rf_distance(other, source, false, &distance)); CHECK(!distance);
  crpx_tree_free(tree); crpx_tree_free(other);
  crpx_tree_node_ids_free(order); crpx_tree_free(source); crpx_distance_matrix_free(distances);
  OK(crpx_distance_matrix_new(3, &distances));
  OK(crpx_distance_matrix_set(distances, 0, 1, 1));
  OK(crpx_distance_matrix_set(distances, 0, 2, 1));
  OK(crpx_distance_matrix_set(distances, 1, 2, 5));
  OK(crpx_tree_neighbour_joining(distances, NULL, true, &tree));
  negative = false;
  for (i = 1; i < crpx_tree_node_count(tree); ++i) {
    OK(crpx_tree_branch_length(tree, (crpx_node_id) i, &length));
    if (length < 0) negative = true;
  }
  CHECK(negative);
  crpx_tree_free(tree);
  OK(crpx_tree_neighbour_joining(distances, NULL, false, &tree));
  for (i = 1; i < crpx_tree_node_count(tree); ++i) {
    OK(crpx_tree_branch_length(tree, (crpx_node_id) i, &length)); CHECK(length >= 0);
  }
  crpx_tree_free(tree);
  OK(crpx_distance_matrix_set(distances, 0, 1, NAN));
  CHECK(crpx_tree_upgma(distances, NULL, &tree) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
  crpx_distance_matrix_free(distances);
  OK(crpx_distance_matrix_new(2, &distances));
  OK(crpx_distance_matrix_set(distances, 0, 1, 6));
  OK(crpx_tree_bionj(distances, NULL, true, &tree));
  OK(crpx_tree_patristic_distance(tree, find_label(tree, "0"), find_label(tree, "1"), NAN, &length));
  CHECK(length == 6);
  crpx_tree_free(tree); crpx_distance_matrix_free(distances);
  OK(crpx_distance_matrix_new(1, &distances));
  OK(crpx_tree_upgma(distances, NULL, &tree)); CHECK(crpx_tree_node_count(tree) == 1);
  crpx_tree_free(tree);
  OK(crpx_tree_bionj(distances, NULL, true, &tree)); CHECK(crpx_tree_node_count(tree) == 1);
  crpx_tree_free(tree); crpx_distance_matrix_free(distances);
  /* A three-tip cluster must contribute three times a single tip in UPGMA. */
  OK(crpx_distance_matrix_new(4, &distances));
  OK(crpx_distance_matrix_set(distances, 0, 1, 2));
  OK(crpx_distance_matrix_set(distances, 0, 2, 4));
  OK(crpx_distance_matrix_set(distances, 1, 2, 4));
  OK(crpx_distance_matrix_set(distances, 0, 3, 10));
  OK(crpx_distance_matrix_set(distances, 1, 3, 10));
  OK(crpx_distance_matrix_set(distances, 2, 3, 6));
  OK(crpx_tree_upgma(distances, NULL, &tree));
  OK(crpx_tree_patristic_distance(tree, find_label(tree, "0"), find_label(tree, "3"), NAN, &length));
  CHECK(fabs(length - 26.0 / 3.0) < 1e-12);
  crpx_tree_free(tree);
  OK(crpx_tree_single_linkage(distances, NULL, &tree));
  OK(crpx_tree_patristic_distance(tree, find_label(tree, "0"), find_label(tree, "3"), NAN, &length));
  CHECK(length == 6);
  crpx_tree_free(tree);
  /* Gascuel's variance difference V_bk - V_ak gives lambda=1 for the first AB join.
   * Its fitted AC and AD paths are 4 and 6; reversing that difference gives 3 and 7.
   */
  OK(crpx_distance_matrix_set(distances, 0, 1, 2));
  OK(crpx_distance_matrix_set(distances, 0, 2, 4));
  OK(crpx_distance_matrix_set(distances, 0, 3, 6));
  OK(crpx_distance_matrix_set(distances, 1, 2, 5));
  OK(crpx_distance_matrix_set(distances, 1, 3, 9));
  OK(crpx_distance_matrix_set(distances, 2, 3, 8));
  OK(crpx_tree_bionj(distances, NULL, true, &tree));
  OK(crpx_tree_patristic_distance(tree, find_label(tree, "0"), find_label(tree, "2"), NAN, &length));
  CHECK(length == 4);
  OK(crpx_tree_patristic_distance(tree, find_label(tree, "0"), find_label(tree, "3"), NAN, &length));
  CHECK(length == 6);
  crpx_tree_free(tree); crpx_distance_matrix_free(distances);
}

static void
test_collections(void)
{
  const char *nexus = "[comment] #NEXUS\nBEGIN TAXA; DIMENSIONS NTAX=3; TAXLABELS A B C; END;\n"
    "Begin Trees; Translate 1 'A semicolon;O''Brien', 2 B, 3 C;\n"
    "TREE * first = [&R] (1,(2,3));\nUTREE second = ((3,2),1); END;";
  const char *bad_nexus[] = {
    "#NEXUS BEGIN TREES; TREE x = (A,B);",
    "#NEXUS BEGIN TREES; TRANSLATE 1 A, 1 B; END;",
    "#NEXUS BEGIN TREES; TRANSLATE 1 A, 2 B; TREE x = (1,3); END;",
    "#NEXUS BEGIN TREES; TREE = (A,B); END;",
    "#NEXUS BEGIN TREES; TREE x (A,B); END;"
  };
  crpx_tree_space *space, *copy, *selected, *read;
  crpx_tree *tree;
  size_t offset, distance, i;
  uint64_t weight;
  double frequency, sum;
  const char *path = "test-tree-space-temporary.nwk";
  char *label;
  FILE *file;
  OK(crpx_tree_space_new(&space));
  OK(crpx_tree_space_add_trees(space, nexus, &offset));
  CHECK(offset == strlen(nexus) && crpx_tree_space_size(space) == 2);
  CHECK(find_label(crpx_tree_space_get(space, 0), "A semicolon;O'Brien") != CRPX_NODE_NONE);
  OK(crpx_tree_space_compact(space, true));
  CHECK(crpx_tree_space_size(space) == 1 && crpx_tree_space_total_weight(space) == 2);
  OK(crpx_tree_space_weight(space, 0, &weight)); CHECK(weight == 2);
  OK(crpx_tree_space_frequency(space, 0, &frequency)); CHECK(frequency == 1);
  for (i = 0; i < sizeof(bad_nexus) / sizeof(*bad_nexus); ++i) {
    CHECK(crpx_tree_space_add_trees(space, bad_nexus[i], &offset) == CRPX_ERROR_PARSE);
    CHECK(offset <= strlen(bad_nexus[i]) && crpx_tree_space_size(space) == 1);
  }
  CHECK(crpx_tree_space_add_newick(space, "(A,B); (C,D)", &offset) == CRPX_ERROR_PARSE);
  CHECK(crpx_tree_space_size(space) == 1);
  OK(crpx_tree_space_write_file(space, path));
  OK(crpx_tree_space_read_file(path, &read, &offset));
  CHECK(crpx_tree_space_size(read) == 1 && crpx_tree_space_total_weight(read) == 1);
  OK(crpx_tree_rf_distance(crpx_tree_space_get(read, 0), crpx_tree_space_get(space, 0), true, &distance));
  CHECK(!distance);
  crpx_tree_space_free(read);
  CHECK(remove(path) == 0);
  crpx_tree_space_free(space);
  OK(crpx_tree_space_new(&space));
  OK(crpx_tree_space_add_trees(space, "#taxon; '#NEXUS';", NULL));
  CHECK(crpx_tree_space_size(space) == 2);
  CHECK(!strcmp(crpx_tree_label(crpx_tree_space_get(space, 0), 0), "#taxon"));
  crpx_tree_space_free(space);
  OK(crpx_tree_space_new(&space));
  OK(crpx_tree_space_add_newick(space, "[;] (A,(B,C)); ((C,B),A); (B,(A,C)); (C,(A,B)); (A,(B,C)); [end]", NULL));
  CHECK(crpx_tree_space_size(space) == 5);
  OK(crpx_tree_space_clone(space, &copy));
  OK(crpx_tree_space_select(space, 1, 2, &selected));
  CHECK(crpx_tree_space_size(selected) == 2);
  crpx_tree_space_free(selected);
  OK(crpx_tree_space_merge(space, space, 2));
  CHECK(crpx_tree_space_size(space) == 10 && crpx_tree_space_total_weight(space) == 15);
  OK(crpx_tree_space_compact(space, true));
  CHECK(crpx_tree_space_size(space) == 3);
  sum = 0;
  for (i = 0; i < 3; ++i) { OK(crpx_tree_space_frequency(space, i, &frequency)); sum += frequency; }
  CHECK(fabs(sum - 1) < 1e-12);
  crpx_tree_space_free(space); crpx_tree_space_free(copy);
  OK(crpx_tree_space_new(&space));
  OK(crpx_tree_from_newick("A;", &tree, NULL));
  OK(crpx_tree_space_add(space, tree, UINT64_MAX));
  CHECK(crpx_tree_space_add(space, tree, 1) == CRPX_ERROR_OVERFLOW);
  CHECK(crpx_tree_space_merge(space, space, 1) == CRPX_ERROR_OVERFLOW);
  CHECK(crpx_tree_space_total_weight(space) == UINT64_MAX);
  crpx_tree_free(tree); crpx_tree_space_free(space);
  /* A quoted semicolon and a doubled quote straddle the stream's 4096-byte blocks. */
  label = (char *) malloc(6001);
  CHECK(label);
  memset(label, 'x', 6000); label[4092] = ';'; label[4093] = '\''; label[6000] = '\0';
  OK(crpx_tree_new(&tree));
  { crpx_node_id root; OK(crpx_tree_add_node(tree, CRPX_NODE_NONE, &root)); OK(crpx_tree_set_label(tree, root, label)); }
  OK(crpx_tree_space_new(&space));
  OK(crpx_tree_space_add(space, tree, 1));
  OK(crpx_tree_space_write_file(space, path));
  OK(crpx_tree_space_read_file(path, &read, &offset));
  CHECK(!strcmp(crpx_tree_label(crpx_tree_space_get(read, 0), 0), label));
  crpx_tree_space_free(read); crpx_tree_space_free(space); crpx_tree_free(tree);
  CHECK(remove(path) == 0);
  free(label);
  file = fopen(path, "wb"); CHECK(file);
  CHECK(fputs(nexus, file) >= 0 && fclose(file) == 0);
  OK(crpx_tree_space_read_file(path, &read, &offset));
  CHECK(crpx_tree_space_size(read) == 2);
  crpx_tree_space_free(read);
  CHECK(remove(path) == 0);
  file = fopen(path, "wb"); CHECK(file);
  CHECK(fputs("(A,B); (C,D)", file) >= 0 && fclose(file) == 0);
  CHECK(crpx_tree_space_read_file(path, &read, &offset) == CRPX_ERROR_PARSE && !read);
  CHECK(remove(path) == 0);
}

int
main(void)
{
  test_rf();
  test_distances();
  test_spr_exhaustive();
  test_randomisation();
  test_distance_construction();
  test_collections();
  return 0;
}
