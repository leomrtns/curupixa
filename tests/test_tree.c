/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 */
#include "tree.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(test) do { if (!(test)) { \
  fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #test); abort(); \
} } while (0)
#define OK(expression) CHECK((expression) == CRPX_OK)

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
check_length(const crpx_tree *tree, crpx_node_id node, double expected)
{
  double length;
  OK(crpx_tree_branch_length(tree, node, &length));
  CHECK(fabs(length - expected) <= 4.0 * DBL_EPSILON * fabs(expected));
}

static void
test_construction(void)
{
  crpx_tree *tree, *clone;
  crpx_node_id root, left, right, child, untouched = 17;
  uint64_t revision;
  double length = 19.0;
  char *text = (char *) 1;
  OK(crpx_tree_new(&tree));
  CHECK(crpx_tree_node_count(tree) == 0 && crpx_tree_root(tree) == CRPX_NODE_NONE);
  OK(crpx_tree_validate(tree));
  CHECK(crpx_tree_to_newick(tree, &text) == CRPX_ERROR_TOPOLOGY && !text);
  CHECK(crpx_tree_add_node(tree, 0, &untouched) == CRPX_ERROR_OUT_OF_RANGE && untouched == 17);
  OK(crpx_tree_reserve(tree, 100));
  CHECK(crpx_tree_storage_size(tree) <= 100 * 12 + 128);
  OK(crpx_tree_add_node(tree, CRPX_NODE_NONE, &root));
  OK(crpx_tree_add_node(tree, root, &left));
  OK(crpx_tree_add_node(tree, root, &right));
  OK(crpx_tree_add_node(tree, left, &child));
  CHECK(root == 0 && left == 1 && right == 2 && child == 3);
  CHECK(crpx_tree_first_child(tree, root) == left);
  CHECK(crpx_tree_next_sibling(tree, left) == right);
  CHECK(crpx_tree_parent(tree, child) == left);
  CHECK(crpx_tree_leaf_count(tree) == 2);
  CHECK(crpx_tree_parent(tree, 100) == CRPX_NODE_NONE);
  CHECK(crpx_tree_branch_length(tree, left, &length) == CRPX_ERROR_NOT_FOUND && length == 19.0);
  OK(crpx_tree_set_label(tree, child, "A label"));
  OK(crpx_tree_set_branch_length(tree, child, 0.125));
  OK(crpx_tree_set_label(tree, right, ""));
  revision = crpx_tree_revision(tree);
  OK(crpx_tree_set_label(tree, right, ""));
  OK(crpx_tree_set_branch_length(tree, child, 0.125));
  CHECK(crpx_tree_revision(tree) == revision);
  CHECK(crpx_tree_set_branch_length(tree, child, INFINITY) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(crpx_tree_set_branch_length(tree, child, NAN) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(crpx_tree_reserve(tree, SIZE_MAX) == CRPX_ERROR_OVERFLOW);
  OK(crpx_tree_validate(tree));
  OK(crpx_tree_to_newick(tree, &text));
  CHECK(!strcmp(text, "(('A label':0.125),'');"));
  crpx_tree_string_free(text);
  OK(crpx_tree_clone(tree, &clone));
  OK(crpx_tree_set_label(tree, child, "changed"));
  CHECK(!strcmp(crpx_tree_label(clone, child), "A label"));
  OK(crpx_tree_clear_branch_length(tree, child));
  CHECK(crpx_tree_branch_length(tree, child, &length) == CRPX_ERROR_NOT_FOUND);
  check_length(clone, child, 0.125);
  OK(crpx_tree_set_label(tree, child, NULL));
  CHECK(!crpx_tree_label(tree, child));
  crpx_tree_free(tree);
  OK(crpx_tree_validate(clone));
  crpx_tree_free(clone);
  crpx_tree_free(NULL);
  crpx_tree_string_free(NULL);
}

static void
test_newick(void)
{
  const char *examples[] = {
    "((A:1,B:2)X:3,C:4)R:5;",
    "[before] ('a b':1.5e-2,[outer[inner]]'O''Brien':-2.125,C_d:0) 'inner label' [after]:3;",
    "(((),''),,A)root;",
    "A:1e-300;", "A:1.7976931348623157e308;", "A:4.9406564584124654e-324;",
    "A:-0;", ";", "('αβγ',中文);", "(A:0000123.4500,B:.5,C:3.);"
  };
  const char *bad[] = {
    "", " ", "A", "(A,B)", "(A,B;", "A;;", "A;B;", "(A B,C);", "'abc;",
    "A:[missing;", "A[unclosed;", "(A,B))", "A:nan;", "A:inf;", "A:1e99999;",
    "A:1e;", "A:;", "A:+;", "A:1.2.3;", "A:0x10;", "A:1_2;", "(A;B);", "A'B;"
  };
  size_t i, j, offset;
  crpx_tree *tree, *roundtrip;
  char *text, *again;
  for (i = 0; i < sizeof(examples) / sizeof(*examples); ++i) {
    OK(crpx_tree_from_newick(examples[i], &tree, &offset));
    CHECK(offset == strlen(examples[i]));
    OK(crpx_tree_validate(tree));
    OK(crpx_tree_to_newick(tree, &text));
    OK(crpx_tree_from_newick(text, &roundtrip, &offset));
    CHECK(crpx_tree_node_count(tree) == crpx_tree_node_count(roundtrip));
    CHECK(crpx_tree_leaf_count(tree) == crpx_tree_leaf_count(roundtrip));
    for (j = 0; j < crpx_tree_node_count(tree); ++j) {
      const char *a = crpx_tree_label(tree, (crpx_node_id) j);
      const char *b = crpx_tree_label(roundtrip, (crpx_node_id) j);
      double x = 0, y = 0;
      crpx_status sx = crpx_tree_branch_length(tree, (crpx_node_id) j, &x);
      crpx_status sy = crpx_tree_branch_length(roundtrip, (crpx_node_id) j, &y);
      CHECK((!a && !b) || (a && b && !strcmp(a, b)));
      CHECK(sx == sy);
      if (sx == CRPX_OK) CHECK(x == y || fabs(x - y) <= 4 * DBL_EPSILON * fabs(x));
      CHECK(crpx_tree_parent(tree, (crpx_node_id) j) == crpx_tree_parent(roundtrip, (crpx_node_id) j));
    }
    OK(crpx_tree_to_newick(roundtrip, &again));
    CHECK(!strcmp(text, again));
    crpx_tree_string_free(again);
    crpx_tree_string_free(text);
    crpx_tree_free(roundtrip);
    crpx_tree_free(tree);
  }
  for (i = 0; i < sizeof(bad) / sizeof(*bad); ++i) {
    tree = (crpx_tree *) 1;
    CHECK(crpx_tree_from_newick(bad[i], &tree, &offset) == CRPX_ERROR_PARSE);
    CHECK(!tree && offset <= strlen(bad[i]));
  }
  CHECK(crpx_tree_from_newick(NULL, &tree, NULL) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
  CHECK(crpx_tree_new(NULL) == CRPX_ERROR_INVALID_ARGUMENT);
  OK(crpx_tree_from_newick(examples[1], &tree, NULL));
  CHECK(find_label(tree, "a b") != CRPX_NODE_NONE);
  CHECK(find_label(tree, "O'Brien") != CRPX_NODE_NONE);
  CHECK(find_label(tree, "C_d") != CRPX_NODE_NONE);
  check_length(tree, find_label(tree, "inner label"), 3);
  crpx_tree_free(tree);
}

static void
test_mutation(void)
{
  crpx_tree *tree;
  crpx_node_id x, a, b, c, old_root;
  uint64_t revision;
  OK(crpx_tree_from_newick("((A:1,B:2)X:3,C:4)R:5;", &tree, NULL));
  x = find_label(tree, "X"); a = find_label(tree, "A");
  b = find_label(tree, "B"); c = find_label(tree, "C");
  old_root = crpx_tree_root(tree);
  revision = crpx_tree_revision(tree);
  CHECK(crpx_tree_reparent(tree, x, a) == CRPX_ERROR_TOPOLOGY);
  CHECK(crpx_tree_reparent(tree, old_root, x) == CRPX_ERROR_TOPOLOGY);
  CHECK(crpx_tree_reparent(tree, x, x) == CRPX_ERROR_TOPOLOGY);
  CHECK(crpx_tree_reroot(tree, a) == CRPX_ERROR_TOPOLOGY);
  CHECK(crpx_tree_revision(tree) == revision);
  OK(crpx_tree_reroot(tree, x));
  OK(crpx_tree_validate(tree));
  CHECK(crpx_tree_root(tree) == x);
  CHECK(crpx_tree_parent(tree, old_root) == x);
  CHECK(crpx_tree_parent(tree, a) == x);
  CHECK(crpx_tree_parent(tree, b) == x);
  check_length(tree, x, 5);
  check_length(tree, old_root, 3);
  check_length(tree, a, 1);
  check_length(tree, c, 4);
  CHECK(crpx_tree_leaf_count(tree) == 3);
  OK(crpx_tree_reroot(tree, old_root));
  OK(crpx_tree_validate(tree));
  check_length(tree, old_root, 5);
  check_length(tree, x, 3);
  OK(crpx_tree_reparent(tree, a, old_root));
  CHECK(crpx_tree_parent(tree, a) == old_root);
  check_length(tree, a, 1);
  OK(crpx_tree_validate(tree));
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("((A,B)X)R;", &tree, NULL));
  revision = crpx_tree_revision(tree);
  CHECK(crpx_tree_reroot(tree, find_label(tree, "X")) == CRPX_ERROR_TOPOLOGY);
  CHECK(crpx_tree_revision(tree) == revision && crpx_tree_leaf_count(tree) == 2);
  OK(crpx_tree_validate(tree));
  crpx_tree_free(tree);
}

static void
test_workspace(void)
{
  crpx_tree *tree, *clone;
  crpx_tree_workspace *workspace, *other;
  const crpx_node_id *preorder, *postorder, *sisters;
  const uint32_t *depth;
  size_t i, n;
  crpx_node_id a, b, x, root;
  uint8_t flags;
  OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_ALL, &workspace));
  OK(crpx_tree_from_newick("((A,B)X,C,D)R;", &tree, NULL));
  OK(crpx_tree_workspace_refresh(workspace, tree));
  CHECK(crpx_tree_workspace_is_current(workspace, tree));
  n = crpx_tree_node_count(tree);
  CHECK(crpx_tree_workspace_node_count(workspace) == n);
  preorder = crpx_tree_workspace_preorder(workspace);
  postorder = crpx_tree_workspace_postorder(workspace);
  sisters = crpx_tree_workspace_sisters(workspace);
  depth = crpx_tree_workspace_depth(workspace);
  root = crpx_tree_root(tree);
  a = find_label(tree, "A"); b = find_label(tree, "B"); x = find_label(tree, "X");
  CHECK(preorder[0] == root && postorder[n - 1] == root);
  CHECK(sisters[a] == b && sisters[b] == a && sisters[x] == CRPX_NODE_NONE);
  CHECK(depth[root] == 0 && depth[x] == 1 && depth[a] == 2);
  for (i = 0; i < n; ++i) CHECK(preorder[i] == i);
  OK(crpx_tree_workspace_set_updated(workspace, tree, a, CRPX_TREE_UPDATED_BOTH));
  OK(crpx_tree_workspace_refresh(workspace, tree));
  OK(crpx_tree_workspace_get_updated(workspace, tree, a, &flags));
  CHECK(flags == CRPX_TREE_UPDATED_BOTH);
  OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_ALL, &other));
  OK(crpx_tree_workspace_refresh(other, tree));
  OK(crpx_tree_workspace_get_updated(other, tree, a, &flags));
  CHECK(!flags);
  OK(crpx_tree_set_branch_length(tree, a, 0.3));
  CHECK(!crpx_tree_workspace_is_current(workspace, tree));
  CHECK(crpx_tree_workspace_get_updated(workspace, tree, a, &flags) == CRPX_ERROR_TOPOLOGY);
  OK(crpx_tree_workspace_refresh(workspace, tree));
  OK(crpx_tree_workspace_get_updated(workspace, tree, a, &flags));
  CHECK(!flags);
  CHECK(!crpx_tree_workspace_is_current(other, tree));
  OK(crpx_tree_clone(tree, &clone));
  CHECK(!crpx_tree_workspace_is_current(workspace, clone));
  OK(crpx_tree_workspace_refresh(workspace, clone));
  CHECK(crpx_tree_workspace_is_current(workspace, clone));
  crpx_tree_workspace_reset(workspace);
  CHECK(!crpx_tree_workspace_is_current(workspace, clone));
  CHECK(crpx_tree_workspace_node_count(workspace) == 0);
  OK(crpx_tree_workspace_refresh(workspace, tree));
  CHECK(crpx_tree_workspace_is_current(workspace, tree));
  crpx_tree_workspace_free(workspace);
  crpx_tree_workspace_free(other);
  OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_DEPTH, &workspace));
  OK(crpx_tree_workspace_refresh(workspace, tree));
  CHECK(!crpx_tree_workspace_preorder(workspace) && !crpx_tree_workspace_sisters(workspace));
  CHECK(crpx_tree_workspace_get_updated(workspace, tree, a, &flags) == CRPX_ERROR_UNSUPPORTED);
  crpx_tree_workspace_free(workspace);
  crpx_tree_free(tree);
  crpx_tree_free(clone);
  OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_ALL, &workspace));
  OK(crpx_tree_new(&tree));
  OK(crpx_tree_workspace_refresh(workspace, tree));
  CHECK(crpx_tree_workspace_node_count(workspace) == 0);
  crpx_tree_workspace_free(workspace);
  crpx_tree_free(tree);
}

static void
test_large_trees(void)
{
  const size_t depth = 50000, leaves = 20000;
  char *input, *text;
  size_t i;
  crpx_tree *tree;
  crpx_tree_workspace *workspace;
  input = (char *) malloc(depth * 2 + 3);
  CHECK(input);
  for (i = 0; i < depth; ++i) input[i] = '(';
  input[depth] = 'A';
  for (i = 0; i < depth; ++i) input[depth + 1 + i] = ')';
  input[depth * 2 + 1] = ';';
  input[depth * 2 + 2] = '\0';
  OK(crpx_tree_from_newick(input, &tree, NULL));
  CHECK(crpx_tree_node_count(tree) == depth + 1);
  OK(crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_ALL, &workspace));
  OK(crpx_tree_workspace_refresh(workspace, tree));
  CHECK(crpx_tree_workspace_depth(workspace)[depth] == depth);
  OK(crpx_tree_to_newick(tree, &text));
  CHECK(!strcmp(input, text));
  crpx_tree_string_free(text);
  crpx_tree_workspace_free(workspace);
  crpx_tree_free(tree);
  free(input);
  input = (char *) malloc(leaves * 2 + 3);
  CHECK(input);
  input[0] = '(';
  for (i = 0; i < leaves; ++i) {
    input[1 + i * 2] = 'A';
    input[2 + i * 2] = i + 1 == leaves ? ')' : ',';
  }
  input[leaves * 2 + 1] = ';';
  input[leaves * 2 + 2] = '\0';
  OK(crpx_tree_from_newick(input, &tree, NULL));
  CHECK(crpx_tree_leaf_count(tree) == leaves);
  OK(crpx_tree_to_newick(tree, &text));
  CHECK(!strcmp(input, text));
  crpx_tree_string_free(text);
  crpx_tree_free(tree);
  free(input);
}

int
main(void)
{
  test_construction();
  test_newick();
  test_mutation();
  test_workspace();
  test_large_trees();
  return 0;
}
