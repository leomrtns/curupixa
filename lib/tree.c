/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * General-purpose tree operations reimplemented from biomcmc-lib's topology and Newick modules.
 */
#include "tree.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  crpx_node_id parent, first_child, next_sibling;
} crpx_tree_node;

struct crpx_tree {
  crpx_tree_node *nodes;
  double *lengths;
  char **labels;
  size_t count, capacity;
  crpx_node_id root;
  uint64_t revision, topology_revision;
};

struct crpx_tree_workspace {
  const crpx_tree *tree;
  uint64_t revision;
  uint32_t features;
  size_t count, capacity;
  crpx_node_id *preorder, *postorder, *sisters;
  uint32_t *depth;
  uint8_t *updated;
};

struct crpx_tree_builder {
  const crpx_tree *tree;
  uint64_t revision;
  crpx_node_id *tails;
  size_t capacity;
};

typedef struct {
  char *data;
  size_t size, capacity;
} tree_string;

static bool
valid_node(const crpx_tree *tree, crpx_node_id node)
{
  return tree && node < tree->count;
}

static crpx_status
string_append(tree_string *string, const char *text, size_t length)
{
  size_t needed, capacity;
  char *replacement;
  if (length > SIZE_MAX - string->size - 1) return CRPX_ERROR_OVERFLOW;
  needed = string->size + length + 1;
  if (needed > string->capacity) {
    capacity = string->capacity ? string->capacity : 64;
    while (capacity < needed) {
      if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
      capacity *= 2;
    }
    replacement = (char *) realloc(string->data, capacity);
    if (!replacement) return CRPX_ERROR_ALLOCATION;
    string->data = replacement;
    string->capacity = capacity;
  }
  if (length) memcpy(string->data + string->size, text, length);
  string->size += length;
  string->data[string->size] = '\0';
  return CRPX_OK;
}

crpx_status
crpx_tree_new(crpx_tree **result)
{
  crpx_tree *tree;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  tree = (crpx_tree *) calloc(1, sizeof(*tree));
  if (!tree) return CRPX_ERROR_ALLOCATION;
  tree->root = CRPX_NODE_NONE;
  tree->revision = 1;
  tree->topology_revision = 1;
  *result = tree;
  return CRPX_OK;
}

void
crpx_tree_free(crpx_tree *tree)
{
  size_t i;
  if (!tree) return;
  if (tree->labels) for (i = 0; i < tree->count; ++i) free(tree->labels[i]);
  free(tree->labels);
  free(tree->lengths);
  free(tree->nodes);
  free(tree);
}

crpx_status
crpx_tree_reserve(crpx_tree *tree, size_t capacity)
{
  crpx_tree_node *nodes;
  double *lengths = NULL;
  char **labels = NULL;
  size_t i;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity <= tree->capacity) return CRPX_OK;
  if (capacity > UINT32_MAX || capacity > SIZE_MAX / sizeof(*nodes) ||
      capacity > SIZE_MAX / sizeof(*lengths) || capacity > SIZE_MAX / sizeof(*labels))
    return CRPX_ERROR_OVERFLOW;
  nodes = (crpx_tree_node *) malloc(capacity * sizeof(*nodes));
  if (!nodes) return CRPX_ERROR_ALLOCATION;
  if (tree->lengths) {
    lengths = (double *) malloc(capacity * sizeof(*lengths));
    if (!lengths) { free(nodes); return CRPX_ERROR_ALLOCATION; }
    if (tree->count) memcpy(lengths, tree->lengths, tree->count * sizeof(*lengths));
    for (i = tree->count; i < capacity; ++i) lengths[i] = NAN;
  }
  if (tree->labels) {
    labels = (char **) calloc(capacity, sizeof(*labels));
    if (!labels) { free(lengths); free(nodes); return CRPX_ERROR_ALLOCATION; }
    if (tree->count) memcpy(labels, tree->labels, tree->count * sizeof(*labels));
  }
  if (tree->count) memcpy(nodes, tree->nodes, tree->count * sizeof(*nodes));
  free(tree->nodes);
  free(tree->lengths);
  free(tree->labels);
  tree->nodes = nodes;
  tree->lengths = lengths;
  tree->labels = labels;
  tree->capacity = capacity;
  return CRPX_OK;
}

static crpx_status
add_node_at_tail(crpx_tree *tree, crpx_node_id parent, crpx_node_id tail, crpx_node_id *result)
{
  crpx_status status;
  size_t capacity;
  crpx_node_id node;
  if (!tree || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  if ((!tree->count && parent != CRPX_NODE_NONE) ||
      (tree->count && !valid_node(tree, parent))) return CRPX_ERROR_OUT_OF_RANGE;
  if (tree->count == UINT32_MAX || tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  if (tree->count == tree->capacity) {
    capacity = tree->capacity ? tree->capacity : 16;
    if (capacity <= tree->count) capacity = capacity > UINT32_MAX / 2 ? UINT32_MAX : capacity * 2;
    status = crpx_tree_reserve(tree, capacity);
    if (status != CRPX_OK) return status;
  }
  node = (crpx_node_id) tree->count;
  tree->nodes[node].parent = parent;
  tree->nodes[node].first_child = CRPX_NODE_NONE;
  tree->nodes[node].next_sibling = CRPX_NODE_NONE;
  if (parent == CRPX_NODE_NONE) tree->root = node;
  else if (tail == CRPX_NODE_NONE) tree->nodes[parent].first_child = node;
  else tree->nodes[tail].next_sibling = node;
  ++tree->count;
  ++tree->revision;
  ++tree->topology_revision;
  *result = node;
  return CRPX_OK;
}

crpx_status
crpx_tree_add_node(crpx_tree *tree, crpx_node_id parent, crpx_node_id *result)
{
  crpx_node_id tail = CRPX_NODE_NONE;
  if (!tree || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  if (tree->count && !valid_node(tree, parent)) return CRPX_ERROR_OUT_OF_RANGE;
  if (valid_node(tree, parent)) {
    tail = tree->nodes[parent].first_child;
    if (tail != CRPX_NODE_NONE)
      while (tree->nodes[tail].next_sibling != CRPX_NODE_NONE) tail = tree->nodes[tail].next_sibling;
  }
  return add_node_at_tail(tree, parent, tail, result);
}

crpx_status crpx_tree_builder_new(crpx_tree_builder **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = (crpx_tree_builder *)calloc(1,sizeof(**result));
  return *result ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void crpx_tree_builder_free(crpx_tree_builder *builder)
{ if (builder) { free(builder->tails); free(builder); } }
void crpx_tree_builder_reset(crpx_tree_builder *builder)
{ if (builder) { builder->tree = NULL; builder->revision = 0; } }
size_t crpx_tree_builder_storage_size(const crpx_tree_builder *builder)
{ return builder ? sizeof(*builder)+builder->capacity*sizeof(crpx_node_id) : 0; }
crpx_status crpx_tree_builder_append(crpx_tree_builder *builder, crpx_tree *tree,
                                    crpx_node_id parent, crpx_node_id *result)
{
  if (!builder || !tree || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  if ((!tree->count && parent != CRPX_NODE_NONE) || (tree->count && !valid_node(tree,parent)))
    return CRPX_ERROR_OUT_OF_RANGE;
  if (tree->count == UINT32_MAX) return CRPX_ERROR_OVERFLOW;
  if (builder->tree != tree || builder->revision != tree->topology_revision) {
    crpx_status status = crpx_tree_validate(tree);
    if (status != CRPX_OK) return status;
    /* Mark stale before reallocating: failed refresh cannot certify old tails. */
    crpx_tree_builder_reset(builder);
  }
  if (tree->count+1 > builder->capacity) {
    size_t capacity = builder->capacity ? builder->capacity : 16;
    while (capacity < tree->count+1) {
      if (capacity > UINT32_MAX/2) { capacity = tree->count+1; break; }
      capacity *= 2;
    }
    if (capacity > SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
    crpx_node_id *tails = (crpx_node_id *)realloc(builder->tails,capacity*sizeof(*tails));
    if (!tails) return CRPX_ERROR_ALLOCATION;
    builder->tails = tails; builder->capacity = capacity;
  }
  if (!builder->tree) {
    for (size_t i=0;i<tree->count;++i) {
      crpx_node_id child = tree->nodes[i].first_child;
      if (child != CRPX_NODE_NONE)
        while (tree->nodes[child].next_sibling != CRPX_NODE_NONE) child = tree->nodes[child].next_sibling;
      builder->tails[i] = child;
    }
    builder->tree = tree; builder->revision = tree->topology_revision;
  }
  crpx_node_id added, tail = parent == CRPX_NODE_NONE ? CRPX_NODE_NONE : builder->tails[parent];
  crpx_status status = add_node_at_tail(tree,parent,tail,&added);
  if (status != CRPX_OK) return status;
  builder->tails[added] = CRPX_NODE_NONE;
  if (parent != CRPX_NODE_NONE) builder->tails[parent] = added;
  builder->revision = tree->topology_revision; *result = added;
  return CRPX_OK;
}

size_t crpx_tree_node_count(const crpx_tree *tree) { return tree ? tree->count : 0; }
crpx_node_id crpx_tree_root(const crpx_tree *tree) { return tree ? tree->root : CRPX_NODE_NONE; }
uint64_t crpx_tree_revision(const crpx_tree *tree) { return tree ? tree->revision : 0; }
crpx_node_id crpx_tree_parent(const crpx_tree *tree, crpx_node_id node)
{ return valid_node(tree, node) ? tree->nodes[node].parent : CRPX_NODE_NONE; }
crpx_node_id crpx_tree_first_child(const crpx_tree *tree, crpx_node_id node)
{ return valid_node(tree, node) ? tree->nodes[node].first_child : CRPX_NODE_NONE; }
crpx_node_id crpx_tree_next_sibling(const crpx_tree *tree, crpx_node_id node)
{ return valid_node(tree, node) ? tree->nodes[node].next_sibling : CRPX_NODE_NONE; }

size_t
crpx_tree_leaf_count(const crpx_tree *tree)
{
  size_t i, count = 0;
  if (tree) for (i = 0; i < tree->count; ++i) count += tree->nodes[i].first_child == CRPX_NODE_NONE;
  return count;
}

size_t
crpx_tree_storage_size(const crpx_tree *tree)
{
  size_t i, size;
  if (!tree) return 0;
  size = sizeof(*tree) + tree->capacity * sizeof(*tree->nodes);
  if (tree->lengths) size += tree->capacity * sizeof(*tree->lengths);
  if (tree->labels) {
    size += tree->capacity * sizeof(*tree->labels);
    for (i = 0; i < tree->count; ++i) if (tree->labels[i]) size += strlen(tree->labels[i]) + 1;
  }
  return size;
}

crpx_status
crpx_tree_validate(const crpx_tree *tree)
{
  uint8_t *seen;
  crpx_node_id *stack, node, child;
  size_t pending = 0, visited = 0, siblings;
  crpx_status status = CRPX_OK;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!tree->count) return tree->root == CRPX_NODE_NONE ? CRPX_OK : CRPX_ERROR_TOPOLOGY;
  if (!valid_node(tree, tree->root) || tree->nodes[tree->root].parent != CRPX_NODE_NONE ||
      tree->nodes[tree->root].next_sibling != CRPX_NODE_NONE) return CRPX_ERROR_TOPOLOGY;
  seen = (uint8_t *) calloc(tree->count, sizeof(*seen));
  stack = (crpx_node_id *) malloc(tree->count * sizeof(*stack));
  if (!seen || !stack) { free(seen); free(stack); return CRPX_ERROR_ALLOCATION; }
  stack[pending++] = tree->root;
  while (pending) {
    node = stack[--pending];
    if (!valid_node(tree, node) || seen[node]) { status = CRPX_ERROR_TOPOLOGY; break; }
    seen[node] = 1;
    ++visited;
    siblings = 0;
    for (child = tree->nodes[node].first_child; child != CRPX_NODE_NONE;
         child = tree->nodes[child].next_sibling) {
      if (!valid_node(tree, child) || tree->nodes[child].parent != node ||
          pending == tree->count || ++siblings > tree->count) { status = CRPX_ERROR_TOPOLOGY; break; }
      stack[pending++] = child;
    }
    if (status != CRPX_OK) break;
  }
  if (visited != tree->count) status = CRPX_ERROR_TOPOLOGY;
  free(stack);
  free(seen);
  return status;
}

static int compare_label_pointers(const void *a, const void *b)
{ return strcmp(*(const char *const *)a,*(const char *const *)b); }
crpx_status crpx_tree_validate_policy(const crpx_tree *tree, uint32_t requirements)
{
  if (requirements & ~UINT32_C(127)) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  const char **names = NULL; size_t named = 0;
  if (requirements & CRPX_TREE_REQUIRE_UNIQUE_TIPS) {
    if (tree->count > SIZE_MAX/sizeof(*names)) return CRPX_ERROR_OVERFLOW;
    names = tree->count ? (const char **)malloc(tree->count*sizeof(*names)) : NULL;
    if (tree->count && !names) return CRPX_ERROR_ALLOCATION;
  }
  for (size_t i=0;i<tree->count;++i) {
    crpx_node_id first = tree->nodes[i].first_child;
    if (first == CRPX_NODE_NONE && (requirements & (CRPX_TREE_REQUIRE_NAMED_TIPS|CRPX_TREE_REQUIRE_UNIQUE_TIPS))) {
      const char *name = tree->labels ? tree->labels[i] : NULL;
      if (!name || !*name) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      if (names) names[named++] = name;
    }
    if (first != CRPX_NODE_NONE && (requirements & CRPX_TREE_REQUIRE_BINARY)) {
      crpx_node_id second = tree->nodes[first].next_sibling;
      if (second == CRPX_NODE_NONE || tree->nodes[second].next_sibling != CRPX_NODE_NONE) {
        status = CRPX_ERROR_UNSUPPORTED; goto done;
      }
    }
    if (i == tree->root && !(requirements & CRPX_TREE_INCLUDE_ROOT_STEM)) continue;
    double length = tree->lengths ? tree->lengths[i] : NAN;
    if (isnan(length)) {
      if (requirements & CRPX_TREE_REQUIRE_COMPLETE_LENGTHS) { status = CRPX_ERROR_NOT_FOUND; goto done; }
    } else if (!isfinite(length) ||
        ((requirements & CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS) && length < 0) ||
        ((requirements & CRPX_TREE_REQUIRE_POSITIVE_LENGTHS) && length <= 0)) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
  }
  if (named > 1) {
    qsort(names,named,sizeof(*names),compare_label_pointers);
    for (size_t i=1;i<named;++i) if (!strcmp(names[i-1],names[i])) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
  }
done:
  free(names); return status;
}

crpx_status
crpx_tree_set_label(crpx_tree *tree, crpx_node_id node, const char *label)
{
  size_t length;
  char *copy = NULL;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  if ((!label && (!tree->labels || !tree->labels[node])) ||
      (label && tree->labels && tree->labels[node] && !strcmp(label, tree->labels[node]))) return CRPX_OK;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  if (label) {
    length = strlen(label);
    if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
    copy = (char *) malloc(length + 1);
    if (!copy) return CRPX_ERROR_ALLOCATION;
    memcpy(copy, label, length + 1);
    if (!tree->labels) {
      tree->labels = (char **) calloc(tree->capacity, sizeof(*tree->labels));
      if (!tree->labels) { free(copy); return CRPX_ERROR_ALLOCATION; }
    }
  }
  free(tree->labels[node]);
  tree->labels[node] = copy;
  ++tree->revision;
  return CRPX_OK;
}

const char *
crpx_tree_label(const crpx_tree *tree, crpx_node_id node)
{
  return valid_node(tree, node) && tree->labels ? tree->labels[node] : NULL;
}

crpx_status
crpx_tree_set_branch_length(crpx_tree *tree, crpx_node_id node, double length)
{
  size_t i;
  if (!tree || !isfinite(length)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  if (tree->lengths && tree->lengths[node] == length &&
      signbit(tree->lengths[node]) == signbit(length)) return CRPX_OK;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  if (!tree->lengths) {
    tree->lengths = (double *) malloc(tree->capacity * sizeof(*tree->lengths));
    if (!tree->lengths) return CRPX_ERROR_ALLOCATION;
    for (i = 0; i < tree->capacity; ++i) tree->lengths[i] = NAN;
  }
  tree->lengths[node] = length;
  ++tree->revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_clear_branch_length(crpx_tree *tree, crpx_node_id node)
{
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  if (!tree->lengths || isnan(tree->lengths[node])) return CRPX_OK;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  tree->lengths[node] = NAN;
  ++tree->revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_branch_length(const crpx_tree *tree, crpx_node_id node, double *length)
{
  if (!tree || !length) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  if (!tree->lengths || isnan(tree->lengths[node])) return CRPX_ERROR_NOT_FOUND;
  *length = tree->lengths[node];
  return CRPX_OK;
}

crpx_status
crpx_tree_clone(const crpx_tree *tree, crpx_tree **result)
{
  crpx_tree *copy;
  crpx_status status;
  size_t i;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_new(&copy);
  if (status != CRPX_OK) return status;
  status = crpx_tree_reserve(copy, tree->count);
  if (status != CRPX_OK) { crpx_tree_free(copy); return status; }
  copy->count = tree->count;
  copy->root = tree->root;
  if (tree->count) memcpy(copy->nodes, tree->nodes, tree->count * sizeof(*tree->nodes));
  for (i = 0; i < tree->count; ++i) {
    if (tree->labels && tree->labels[i]) {
      status = crpx_tree_set_label(copy, (crpx_node_id) i, tree->labels[i]);
      if (status != CRPX_OK) { crpx_tree_free(copy); return status; }
    }
    if (tree->lengths && !isnan(tree->lengths[i])) {
      status = crpx_tree_set_branch_length(copy, (crpx_node_id) i, tree->lengths[i]);
      if (status != CRPX_OK) { crpx_tree_free(copy); return status; }
    }
  }
  *result = copy;
  return CRPX_OK;
}

static void
detach_child(crpx_tree *tree, crpx_node_id parent, crpx_node_id node)
{
  crpx_node_id child = tree->nodes[parent].first_child;
  if (child == node) tree->nodes[parent].first_child = tree->nodes[node].next_sibling;
  else {
    while (tree->nodes[child].next_sibling != node) child = tree->nodes[child].next_sibling;
    tree->nodes[child].next_sibling = tree->nodes[node].next_sibling;
  }
  tree->nodes[node].next_sibling = CRPX_NODE_NONE;
}

static void
append_child(crpx_tree *tree, crpx_node_id parent, crpx_node_id node)
{
  crpx_node_id child = tree->nodes[parent].first_child;
  if (child == CRPX_NODE_NONE) tree->nodes[parent].first_child = node;
  else {
    while (tree->nodes[child].next_sibling != CRPX_NODE_NONE) child = tree->nodes[child].next_sibling;
    tree->nodes[child].next_sibling = node;
  }
  tree->nodes[node].parent = parent;
  tree->nodes[node].next_sibling = CRPX_NODE_NONE;
}

static void
replace_child(crpx_tree *tree, crpx_node_id parent, crpx_node_id old_node, crpx_node_id new_node)
{
  crpx_node_id previous = tree->nodes[parent].first_child;
  if (previous == old_node) tree->nodes[parent].first_child = new_node;
  else {
    while (tree->nodes[previous].next_sibling != old_node) previous = tree->nodes[previous].next_sibling;
    tree->nodes[previous].next_sibling = new_node;
  }
}

crpx_status
crpx_tree_swap_subtrees(crpx_tree *tree, crpx_node_id first, crpx_node_id second)
{
  crpx_node_id first_parent, second_parent, first_next, second_next, ancestor;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, first) || !valid_node(tree, second)) return CRPX_ERROR_OUT_OF_RANGE;
  if (first == tree->root || second == tree->root) return CRPX_ERROR_TOPOLOGY;
  if (first == second) return CRPX_OK;
  first_parent = tree->nodes[first].parent;
  second_parent = tree->nodes[second].parent;
  if (first_parent == second_parent) return CRPX_OK;
  /* The uncle/nephew relationship used by NNI proves disjointness in constant time. */
  if (!((tree->nodes[first_parent].parent == second_parent && first_parent != second) ||
        (tree->nodes[second_parent].parent == first_parent && second_parent != first))) {
    for (ancestor = first_parent; ancestor != CRPX_NODE_NONE; ancestor = tree->nodes[ancestor].parent)
      if (ancestor == second) return CRPX_ERROR_TOPOLOGY;
    for (ancestor = second_parent; ancestor != CRPX_NODE_NONE; ancestor = tree->nodes[ancestor].parent)
      if (ancestor == first) return CRPX_ERROR_TOPOLOGY;
  }
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  first_next = tree->nodes[first].next_sibling;
  second_next = tree->nodes[second].next_sibling;
  replace_child(tree, first_parent, first, second);
  replace_child(tree, second_parent, second, first);
  tree->nodes[first].parent = second_parent;
  tree->nodes[first].next_sibling = second_next;
  tree->nodes[second].parent = first_parent;
  tree->nodes[second].next_sibling = first_next;
  ++tree->revision;
  ++tree->topology_revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_spr(crpx_tree *tree, crpx_node_id prune, crpx_node_id regraft, double graft_fraction)
{
  crpx_node_id parent, grandparent, sister, first, second, ancestor, target_parent, target_next;
  double collapsed_length = NAN, target_length = NAN, stem_length = NAN;
  if (!tree || !isfinite(graft_fraction) || graft_fraction < 0 || graft_fraction > 1)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, prune) || !valid_node(tree, regraft)) return CRPX_ERROR_OUT_OF_RANGE;
  if (prune == tree->root) return CRPX_ERROR_TOPOLOGY;
  parent = tree->nodes[prune].parent;
  first = tree->nodes[parent].first_child;
  second = tree->nodes[first].next_sibling;
  if (second == CRPX_NODE_NONE || tree->nodes[second].next_sibling != CRPX_NODE_NONE)
    return CRPX_ERROR_TOPOLOGY;
  sister = first == prune ? second : first;
  if (regraft == parent || regraft == sister) return CRPX_ERROR_TOPOLOGY;
  for (ancestor = regraft; ancestor != CRPX_NODE_NONE; ancestor = tree->nodes[ancestor].parent)
    if (ancestor == prune) return CRPX_ERROR_TOPOLOGY;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  grandparent = tree->nodes[parent].parent;
  if (tree->lengths) {
    collapsed_length = tree->lengths[parent] + tree->lengths[sister];
    if (grandparent != CRPX_NODE_NONE && isinf(collapsed_length)) return CRPX_ERROR_OVERFLOW;
    target_length = tree->lengths[regraft];
    stem_length = tree->lengths[tree->root];
  }
  if (grandparent == CRPX_NODE_NONE) {
    tree->root = sister;
    tree->nodes[sister].parent = CRPX_NODE_NONE;
    tree->nodes[sister].next_sibling = CRPX_NODE_NONE;
    if (tree->lengths) tree->lengths[sister] = stem_length;
  }
  else {
    replace_child(tree, grandparent, parent, sister);
    tree->nodes[sister].parent = grandparent;
    tree->nodes[sister].next_sibling = tree->nodes[parent].next_sibling;
    if (tree->lengths) tree->lengths[sister] = collapsed_length;
  }
  target_parent = tree->nodes[regraft].parent;
  target_next = tree->nodes[regraft].next_sibling;
  if (target_parent == CRPX_NODE_NONE) {
    tree->root = parent;
    tree->nodes[parent].parent = CRPX_NODE_NONE;
    tree->nodes[parent].next_sibling = CRPX_NODE_NONE;
    if (tree->lengths) { tree->lengths[parent] = stem_length; tree->lengths[regraft] = 0; }
  }
  else {
    replace_child(tree, target_parent, regraft, parent);
    tree->nodes[parent].parent = target_parent;
    tree->nodes[parent].next_sibling = target_next;
    if (tree->lengths) {
      tree->lengths[parent] = target_length * graft_fraction;
      tree->lengths[regraft] = target_length * (1 - graft_fraction);
    }
  }
  tree->nodes[parent].first_child = prune;
  tree->nodes[prune].parent = parent;
  tree->nodes[prune].next_sibling = regraft;
  tree->nodes[regraft].parent = parent;
  tree->nodes[regraft].next_sibling = CRPX_NODE_NONE;
  ++tree->revision;
  ++tree->topology_revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_reparent(crpx_tree *tree, crpx_node_id node, crpx_node_id parent)
{
  crpx_node_id ancestor;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node) || !valid_node(tree, parent)) return CRPX_ERROR_OUT_OF_RANGE;
  if (node == tree->root) return CRPX_ERROR_TOPOLOGY;
  for (ancestor = parent; ancestor != CRPX_NODE_NONE; ancestor = tree->nodes[ancestor].parent)
    if (ancestor == node) return CRPX_ERROR_TOPOLOGY;
  if (tree->nodes[node].parent == parent) return CRPX_OK;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  detach_child(tree, tree->nodes[node].parent, node);
  append_child(tree, parent, node);
  ++tree->revision;
  ++tree->topology_revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_reroot(crpx_tree *tree, crpx_node_id node)
{
  crpx_node_id parent, grandparent, current;
  double edge_length = NAN, next_length, stem_length = NAN;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  if (node == tree->root) return CRPX_OK;
  if (tree->nodes[node].first_child == CRPX_NODE_NONE ||
      tree->nodes[tree->nodes[tree->root].first_child].next_sibling == CRPX_NODE_NONE)
    return CRPX_ERROR_TOPOLOGY;
  if (tree->revision == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
  if (tree->lengths) {
    edge_length = tree->lengths[node];
    stem_length = tree->lengths[tree->root];
  }
  current = node;
  parent = tree->nodes[node].parent;
  detach_child(tree, parent, node);
  tree->nodes[node].parent = CRPX_NODE_NONE;
  while (parent != CRPX_NODE_NONE) {
    grandparent = tree->nodes[parent].parent;
    if (grandparent != CRPX_NODE_NONE) detach_child(tree, grandparent, parent);
    append_child(tree, current, parent);
    if (tree->lengths) {
      next_length = tree->lengths[parent];
      tree->lengths[parent] = edge_length;
      edge_length = next_length;
    }
    current = parent;
    parent = grandparent;
  }
  tree->root = node;
  if (tree->lengths) tree->lengths[node] = stem_length;
  ++tree->revision;
  ++tree->topology_revision;
  return CRPX_OK;
}

crpx_status
crpx_tree_reroot_edge(crpx_tree *tree, crpx_node_id child, double fraction,
                      crpx_tree_edge_root_policy policy, crpx_node_id *result)
{
  if (!tree || !result || !isfinite(fraction) || fraction < 0 || fraction > 1 ||
      (policy != CRPX_TREE_REUSE_BINARY_ROOT && policy != CRPX_TREE_KEEP_OLD_ROOT))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (!valid_node(tree,child)) return CRPX_ERROR_OUT_OF_RANGE;
  if (child == tree->root) return CRPX_ERROR_TOPOLOGY;
  if (tree->revision == UINT64_MAX ||
      (policy == CRPX_TREE_KEEP_OLD_ROOT && tree->count == UINT32_MAX)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  crpx_node_id old_root = tree->root, parent = tree->nodes[child].parent;
  crpx_node_id a = tree->nodes[old_root].first_child, b = tree->nodes[a].next_sibling;
  if (b == CRPX_NODE_NONE) return CRPX_ERROR_UNSUPPORTED;
  if (policy == CRPX_TREE_REUSE_BINARY_ROOT && tree->nodes[b].next_sibling != CRPX_NODE_NONE)
    return CRPX_ERROR_UNSUPPORTED;
  size_t count = tree->count+(policy == CRPX_TREE_KEEP_OLD_ROOT), capacity = tree->capacity;
  if (count > capacity) capacity = capacity > UINT32_MAX/2 ? UINT32_MAX : capacity*2;
  if (capacity < count || capacity > SIZE_MAX/sizeof(crpx_tree_node) ||
      capacity > SIZE_MAX/sizeof(double) || capacity > SIZE_MAX/sizeof(char *)) return CRPX_ERROR_OVERFLOW;
  crpx_tree copy = *tree;
  copy.nodes = (crpx_tree_node *)malloc(capacity*sizeof(*copy.nodes));
  copy.lengths = tree->lengths ? (double *)malloc(capacity*sizeof(*copy.lengths)) : NULL;
  copy.labels = tree->labels ? (char **)calloc(capacity,sizeof(*copy.labels)) : NULL;
  if (!copy.nodes || (tree->lengths && !copy.lengths) || (tree->labels && !copy.labels)) {
    free(copy.nodes); free(copy.lengths); free(copy.labels); return CRPX_ERROR_ALLOCATION;
  }
  memcpy(copy.nodes,tree->nodes,tree->count*sizeof(*copy.nodes));
  if (copy.labels) memcpy(copy.labels,tree->labels,tree->count*sizeof(*copy.labels));
  if (copy.lengths) {
    memcpy(copy.lengths,tree->lengths,tree->count*sizeof(*copy.lengths));
    for (size_t i=tree->count;i<capacity;++i) copy.lengths[i] = NAN;
  }
  copy.capacity = capacity; copy.revision = 1;
  crpx_node_id root = policy == CRPX_TREE_KEEP_OLD_ROOT ? (crpx_node_id)tree->count : old_root;
  double edge = copy.lengths ? copy.lengths[child] : NAN;
  double stem = copy.lengths ? copy.lengths[old_root] : NAN;
  double to_parent = edge*fraction, to_child = edge*(1-fraction);
  if (policy == CRPX_TREE_REUSE_BINARY_ROOT) {
    double joined = copy.lengths ? copy.lengths[a]+copy.lengths[b] : NAN;
    if (isinf(joined)) { status = CRPX_ERROR_OVERFLOW; goto failed; }
    if (parent == old_root) {
      crpx_node_id other = child == a ? b : a;
      to_parent += copy.lengths ? copy.lengths[other] : NAN;
      if (isinf(to_parent)) { status = CRPX_ERROR_OVERFLOW; goto failed; }
      parent = other; a = other; b = child;
    } else if (copy.nodes[a].first_child == CRPX_NODE_NONE) {
      crpx_node_id swap = a; a = b; b = swap;
    }
    copy.root = a; copy.nodes[a].parent = CRPX_NODE_NONE; copy.nodes[a].next_sibling = CRPX_NODE_NONE;
    append_child(&copy,a,b);
    if (copy.lengths) { copy.lengths[a] = stem; copy.lengths[b] = joined; }
    copy.nodes[old_root] = (crpx_tree_node){CRPX_NODE_NONE,CRPX_NODE_NONE,CRPX_NODE_NONE};
  }
  status = crpx_tree_reroot(&copy,parent);
  if (status != CRPX_OK) goto failed;
  detach_child(&copy,parent,child);
  copy.nodes[root] = (crpx_tree_node){CRPX_NODE_NONE,child,CRPX_NODE_NONE};
  copy.nodes[child].parent = root; copy.nodes[child].next_sibling = parent;
  copy.nodes[parent].parent = root; copy.nodes[parent].next_sibling = CRPX_NODE_NONE;
  copy.root = root; copy.count = count;
  if (copy.lengths) { copy.lengths[root] = stem; copy.lengths[parent] = to_parent; copy.lengths[child] = to_child; }
  free(tree->nodes); free(tree->lengths); free(tree->labels);
  tree->nodes = copy.nodes; tree->lengths = copy.lengths; tree->labels = copy.labels;
  tree->count = count; tree->capacity = capacity; tree->root = root;
  ++tree->revision; ++tree->topology_revision;
  *result = root; return CRPX_OK;
failed:
  free(copy.nodes); free(copy.lengths); free(copy.labels); return status;
}

/* Skip whitespace and nested Newick comments, without consulting or changing a process locale. */
static bool
newick_space(unsigned char c)
{
  return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\v' || c == '\f';
}

static crpx_status
skip_newick_space(const char *text, size_t *position)
{
  size_t depth;
  for (;;) {
    while (newick_space((unsigned char) text[*position])) ++*position;
    if (text[*position] != '[') return CRPX_OK;
    depth = 1;
    ++*position;
    while (depth) {
      if (!text[*position]) return CRPX_ERROR_PARSE;
      if (text[*position] == '[') ++depth;
      else if (text[*position] == ']') --depth;
      ++*position;
    }
  }
}

static bool
newick_delimiter(char c)
{
  return c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '[' || c == ']';
}

static crpx_status
parse_label(const char *text, size_t *position, char **result)
{
  size_t begin;
  tree_string label = { NULL, 0, 0 };
  crpx_status status;
  *result = NULL;
  if (text[*position] == '\'') {
    ++*position;
    for (;;) {
      if (!text[*position]) { free(label.data); return CRPX_ERROR_PARSE; }
      if (text[*position] == '\'') {
        ++*position;
        if (text[*position] != '\'') break;
      }
      status = string_append(&label, text + *position, 1);
      if (status != CRPX_OK) { free(label.data); return status; }
      ++*position;
    }
    if (!label.data) {
      status = string_append(&label, "", 0);
      if (status != CRPX_OK) return status;
    }
  }
  else {
    begin = *position;
    while (text[*position] && !newick_delimiter(text[*position]) &&
           !newick_space((unsigned char) text[*position])) {
      if (text[*position] == '\'') return CRPX_ERROR_PARSE;
      ++*position;
    }
    if (begin == *position) return CRPX_OK;
    status = string_append(&label, text + begin, *position - begin);
    if (status != CRPX_OK) { free(label.data); return status; }
  }
  *result = label.data;
  return CRPX_OK;
}

static bool
ascii_digit(char c)
{
  return c >= '0' && c <= '9';
}

/* Decimal conversion deliberately avoids strtod's process-locale dependency and shared errno.
 * Retain 19 significant digits in long double; format/parse round trips are within floating-point rounding.
 */
static crpx_status
parse_length(const char *text, size_t *position, double *result)
{
  bool negative = false, decimal = false, any_digit = false, significant = false, exponent_negative = false;
  unsigned kept = 0;
  long double mantissa = 0.0L, value;
  int64_t fractional = 0, skipped = 0, exponent = 0;
  char c;
  if (text[*position] == '+' || text[*position] == '-') negative = text[(*position)++] == '-';
  for (;;) {
    c = text[*position];
    if (ascii_digit(c)) {
      any_digit = true;
      if (decimal) ++fractional;
      if (c != '0') significant = true;
      if (significant) {
        if (kept < 19) { mantissa = mantissa * 10.0L + (c - '0'); ++kept; }
        else ++skipped;
      }
      ++*position;
    }
    else if (c == '.' && !decimal) { decimal = true; ++*position; }
    else break;
  }
  if (!any_digit) return CRPX_ERROR_PARSE;
  if (text[*position] == 'e' || text[*position] == 'E') {
    ++*position;
    if (text[*position] == '+' || text[*position] == '-') exponent_negative = text[(*position)++] == '-';
    if (!ascii_digit(text[*position])) return CRPX_ERROR_PARSE;
    while (ascii_digit(text[*position])) {
      if (exponent < 1000000) exponent = exponent * 10 + (text[*position] - '0');
      ++*position;
    }
  }
  if (exponent_negative) exponent = -exponent;
  exponent += skipped - fractional;
  value = mantissa == 0.0L ? 0.0L : mantissa * powl(10.0L, (long double) exponent);
  if (negative) value = -value;
  *result = (double) value;
  return isfinite(*result) ? CRPX_OK : CRPX_ERROR_PARSE;
}

static crpx_status
parse_suffix(const char *text, size_t *position, crpx_tree *tree, crpx_node_id node)
{
  char *label;
  double length;
  crpx_status status = skip_newick_space(text, position);
  if (status != CRPX_OK) return status;
  status = parse_label(text, position, &label);
  if (status != CRPX_OK) return status;
  if (label) {
    status = crpx_tree_set_label(tree, node, label);
    free(label);
    if (status != CRPX_OK) return status;
  }
  status = skip_newick_space(text, position);
  if (status != CRPX_OK) return status;
  if (text[*position] == ':') {
    ++*position;
    status = skip_newick_space(text, position);
    if (status != CRPX_OK) return status;
    status = parse_length(text, position, &length);
    if (status != CRPX_OK) return status;
    status = crpx_tree_set_branch_length(tree, node, length);
    if (status != CRPX_OK) return status;
  }
  return skip_newick_space(text, position);
}

typedef struct { crpx_node_id node, tail; } newick_frame;

crpx_status
crpx_tree_from_newick(const char *newick, crpx_tree **result, size_t *error_offset)
{
  crpx_tree *tree = NULL;
  newick_frame *stack = NULL, *replacement;
  size_t position = 0, depth = 0, capacity = 0;
  crpx_node_id parent, node;
  bool expect_subtree = true;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (error_offset) *error_offset = 0;
  if (!newick) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_new(&tree);
  if (status != CRPX_OK) return status;
  for (;;) {
    status = skip_newick_space(newick, &position);
    if (status != CRPX_OK) break;
    if (expect_subtree) {
      if (!newick[position] || (newick[position] == ';' && depth)) { status = CRPX_ERROR_PARSE; break; }
      if (newick[position] == '(' && depth == capacity) {
        size_t next = capacity ? capacity * 2 : 16;
        if (next < capacity || next > SIZE_MAX / sizeof(*stack)) { status = CRPX_ERROR_OVERFLOW; break; }
        replacement = (newick_frame *) realloc(stack, next * sizeof(*stack));
        if (!replacement) { status = CRPX_ERROR_ALLOCATION; break; }
        stack = replacement;
        capacity = next;
      }
      parent = depth ? stack[depth - 1].node : CRPX_NODE_NONE;
      status = add_node_at_tail(tree, parent, depth ? stack[depth - 1].tail : CRPX_NODE_NONE, &node);
      if (status != CRPX_OK) break;
      if (depth) stack[depth - 1].tail = node;
      if (newick[position] == '(') {
        ++position;
        stack[depth].node = node;
        stack[depth++].tail = CRPX_NODE_NONE;
        continue;
      }
      status = parse_suffix(newick, &position, tree, node);
      if (status != CRPX_OK) break;
      expect_subtree = false;
    }
    if (newick[position] == ',' && depth) { ++position; expect_subtree = true; }
    else if (newick[position] == ')' && depth) {
      node = stack[--depth].node;
      ++position;
      status = parse_suffix(newick, &position, tree, node);
      if (status != CRPX_OK) break;
    }
    else if (newick[position] == ';' && !depth) {
      ++position;
      status = skip_newick_space(newick, &position);
      if (status == CRPX_OK && newick[position]) status = CRPX_ERROR_PARSE;
      break;
    }
    else { status = CRPX_ERROR_PARSE; break; }
  }
  free(stack);
  if (error_offset) *error_offset = position;
  if (status != CRPX_OK) { crpx_tree_free(tree); return status; }
  *result = tree;
  return CRPX_OK;
}

static crpx_status
write_suffix(const crpx_tree *tree, crpx_node_id node, tree_string *string)
{
  const char *label = crpx_tree_label(tree, node);
  crpx_status status;
  size_t i;
  bool quote;
  char number[128], decimal[128];
  int count;
  if (label) {
    quote = !*label;
    for (i = 0; label[i]; ++i)
      if (newick_space((unsigned char) label[i]) || newick_delimiter(label[i]) || label[i] == '\'') quote = true;
    if (quote) {
      status = string_append(string, "'", 1);
      if (status != CRPX_OK) return status;
    }
    for (i = 0; label[i]; ++i) {
      status = string_append(string, label + i, 1);
      if (status != CRPX_OK) return status;
      if (quote && label[i] == '\'') {
        status = string_append(string, "'", 1);
        if (status != CRPX_OK) return status;
      }
    }
    if (quote) {
      status = string_append(string, "'", 1);
      if (status != CRPX_OK) return status;
    }
  }
  if (tree->lengths && !isnan(tree->lengths[node])) {
    size_t source = 0, target = 0;
    count = snprintf(number, sizeof(number), "%.17g", tree->lengths[node]);
    if (count < 0 || (size_t) count >= sizeof(number)) return CRPX_ERROR_OVERFLOW;
    /* %g uses ASCII digits and a locale-specific decimal separator. Replace that separator
     * without calling setlocale or localeconv, both of which expose shared locale state.
     */
    while (number[source]) {
      if (ascii_digit(number[source]) || number[source] == '+' || number[source] == '-' ||
          number[source] == 'e' || number[source] == 'E') decimal[target++] = number[source++];
      else {
        decimal[target++] = '.';
        do { ++source; }
        while (number[source] && !ascii_digit(number[source]) && number[source] != 'e' && number[source] != 'E');
      }
    }
    status = string_append(string, ":", 1);
    if (status != CRPX_OK) return status;
    status = string_append(string, decimal, target);
    if (status != CRPX_OK) return status;
  }
  return CRPX_OK;
}

crpx_status
crpx_tree_to_newick(const crpx_tree *tree, char **result)
{
  tree_string string = { NULL, 0, 0 };
  crpx_node_id node;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!tree->count) return CRPX_ERROR_TOPOLOGY;
  status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  node = tree->root;
  for (;;) {
    if (tree->nodes[node].first_child != CRPX_NODE_NONE) {
      status = string_append(&string, "(", 1);
      if (status != CRPX_OK) break;
      node = tree->nodes[node].first_child;
      continue;
    }
    status = write_suffix(tree, node, &string);
    if (status != CRPX_OK) break;
    for (;;) {
      if (node == tree->root) {
        status = string_append(&string, ";", 1);
        goto finished;
      }
      if (tree->nodes[node].next_sibling != CRPX_NODE_NONE) {
        status = string_append(&string, ",", 1);
        if (status != CRPX_OK) goto finished;
        node = tree->nodes[node].next_sibling;
        break;
      }
      node = tree->nodes[node].parent;
      status = string_append(&string, ")", 1);
      if (status != CRPX_OK) goto finished;
      status = write_suffix(tree, node, &string);
      if (status != CRPX_OK) goto finished;
    }
  }
finished:
  if (status != CRPX_OK) { free(string.data); return status; }
  *result = string.data;
  return CRPX_OK;
}

void crpx_tree_string_free(char *string) { free(string); }

crpx_status
crpx_tree_workspace_new(uint32_t features, crpx_tree_workspace **result)
{
  crpx_tree_workspace *workspace;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!features || (features & ~((uint32_t) CRPX_TREE_WORKSPACE_ALL))) return CRPX_ERROR_INVALID_ARGUMENT;
  workspace = (crpx_tree_workspace *) calloc(1, sizeof(*workspace));
  if (!workspace) return CRPX_ERROR_ALLOCATION;
  workspace->features = features;
  *result = workspace;
  return CRPX_OK;
}

void
crpx_tree_workspace_free(crpx_tree_workspace *workspace)
{
  if (!workspace) return;
  free(workspace->preorder);
  free(workspace->postorder);
  free(workspace->depth);
  free(workspace->sisters);
  free(workspace->updated);
  free(workspace);
}

void
crpx_tree_workspace_reset(crpx_tree_workspace *workspace)
{
  if (!workspace) return;
  workspace->tree = NULL;
  workspace->revision = 0;
  workspace->count = 0;
}

bool
crpx_tree_workspace_is_current(const crpx_tree_workspace *workspace, const crpx_tree *tree)
{
  return workspace && tree && workspace->tree == tree && workspace->revision == tree->revision;
}

crpx_status
crpx_tree_workspace_refresh(crpx_tree_workspace *workspace, const crpx_tree *tree)
{
  crpx_node_id *preorder = NULL, *postorder = NULL, *sisters = NULL;
  uint32_t *depths = NULL;
  uint8_t *updated = NULL;
  crpx_node_id node, first, second;
  uint32_t depth = 0;
  size_t pre = 0, post = 0, i;
  crpx_status status;
  if (!workspace || !tree) return CRPX_ERROR_INVALID_ARGUMENT;
  if (crpx_tree_workspace_is_current(workspace, tree)) return CRPX_OK;
  status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  if (tree->count > workspace->capacity) {
    if (workspace->features & CRPX_TREE_WORKSPACE_TRAVERSAL) {
      preorder = (crpx_node_id *) malloc(tree->count * sizeof(*preorder));
      postorder = (crpx_node_id *) malloc(tree->count * sizeof(*postorder));
      if (!preorder || !postorder) goto allocation_failed;
    }
    if (workspace->features & CRPX_TREE_WORKSPACE_DEPTH) {
      depths = (uint32_t *) malloc(tree->count * sizeof(*depths));
      if (!depths) goto allocation_failed;
    }
    if (workspace->features & CRPX_TREE_WORKSPACE_SISTERS) {
      sisters = (crpx_node_id *) malloc(tree->count * sizeof(*sisters));
      if (!sisters) goto allocation_failed;
    }
    if (workspace->features & CRPX_TREE_WORKSPACE_UPDATED) {
      updated = (uint8_t *) malloc(tree->count * sizeof(*updated));
      if (!updated) goto allocation_failed;
    }
    free(workspace->preorder); workspace->preorder = preorder;
    free(workspace->postorder); workspace->postorder = postorder;
    free(workspace->depth); workspace->depth = depths;
    free(workspace->sisters); workspace->sisters = sisters;
    free(workspace->updated); workspace->updated = updated;
    workspace->capacity = tree->count;
  }
  if (workspace->updated && tree->count) memset(workspace->updated, 0, tree->count * sizeof(*workspace->updated));
  if (workspace->sisters) {
    for (i = 0; i < tree->count; ++i) workspace->sisters[i] = CRPX_NODE_NONE;
    for (i = 0; i < tree->count; ++i) {
      first = tree->nodes[i].first_child;
      if (first == CRPX_NODE_NONE) continue;
      second = tree->nodes[first].next_sibling;
      if (second != CRPX_NODE_NONE && tree->nodes[second].next_sibling == CRPX_NODE_NONE) {
        workspace->sisters[first] = second;
        workspace->sisters[second] = first;
      }
    }
  }
  node = tree->root;
  while (node != CRPX_NODE_NONE) {
    if (workspace->preorder) workspace->preorder[pre++] = node;
    if (workspace->depth) workspace->depth[node] = depth;
    if (tree->nodes[node].first_child != CRPX_NODE_NONE) {
      node = tree->nodes[node].first_child;
      ++depth;
      continue;
    }
    for (;;) {
      if (workspace->postorder) workspace->postorder[post++] = node;
      if (tree->nodes[node].next_sibling != CRPX_NODE_NONE) { node = tree->nodes[node].next_sibling; break; }
      node = tree->nodes[node].parent;
      if (node == CRPX_NODE_NONE) break;
      --depth;
    }
  }
  workspace->tree = tree;
  workspace->revision = tree->revision;
  workspace->count = tree->count;
  return CRPX_OK;
allocation_failed:
  free(preorder); free(postorder); free(depths); free(sisters); free(updated);
  return CRPX_ERROR_ALLOCATION;
}

size_t crpx_tree_workspace_node_count(const crpx_tree_workspace *workspace)
{ return workspace ? workspace->count : 0; }
const crpx_node_id *crpx_tree_workspace_preorder(const crpx_tree_workspace *workspace)
{ return workspace ? workspace->preorder : NULL; }
const crpx_node_id *crpx_tree_workspace_postorder(const crpx_tree_workspace *workspace)
{ return workspace ? workspace->postorder : NULL; }
const uint32_t *crpx_tree_workspace_depth(const crpx_tree_workspace *workspace)
{ return workspace ? workspace->depth : NULL; }
const crpx_node_id *crpx_tree_workspace_sisters(const crpx_tree_workspace *workspace)
{ return workspace ? workspace->sisters : NULL; }

crpx_status
crpx_tree_workspace_set_updated(crpx_tree_workspace *workspace, const crpx_tree *tree,
                                crpx_node_id node, uint8_t flags)
{
  if (!workspace || !tree || flags > CRPX_TREE_UPDATED_BOTH) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!(workspace->features & CRPX_TREE_WORKSPACE_UPDATED)) return CRPX_ERROR_UNSUPPORTED;
  if (!crpx_tree_workspace_is_current(workspace, tree)) return CRPX_ERROR_TOPOLOGY;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  workspace->updated[node] = flags;
  return CRPX_OK;
}

crpx_status
crpx_tree_workspace_get_updated(const crpx_tree_workspace *workspace, const crpx_tree *tree,
                                crpx_node_id node, uint8_t *flags)
{
  if (!workspace || !tree || !flags) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!(workspace->features & CRPX_TREE_WORKSPACE_UPDATED)) return CRPX_ERROR_UNSUPPORTED;
  if (!crpx_tree_workspace_is_current(workspace, tree)) return CRPX_ERROR_TOPOLOGY;
  if (!valid_node(tree, node)) return CRPX_ERROR_OUT_OF_RANGE;
  *flags = workspace->updated[node];
  return CRPX_OK;
}
