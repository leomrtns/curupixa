/* Reimplementation of biomcmc-lib splitset distances (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "splitset_distances.h"
#include "hashtable.h"
#include "lowlevel_algorithms.h"
#include <stdlib.h>

struct crpx_splitset { size_t taxa, count; crpx_bipartition_t *splits; };

void crpx_splitset_free(crpx_splitset *s)
{
  if (!s) return;
  for (size_t i = 0; i < s->count; ++i) crpx_bipartition_free(s->splits[i]);
  free(s->splits); free(s);
}
size_t crpx_splitset_count(const crpx_splitset *s) { return s ? s->count : 0; }
const struct crpx_bipartition_struct *crpx_splitset_get(const crpx_splitset *s, size_t i)
{ return s && i < s->count ? s->splits[i] : NULL; }

crpx_status crpx_splitset_from_tree(const crpx_tree *tree, const crpx_node_id *tips, size_t taxa, crpx_splitset **out)
{
  crpx_tree_workspace *w = NULL;
  crpx_hashtable_t unique = NULL;
  crpx_splitset *s = NULL;
  crpx_bipartition_t *descendants = NULL;
  crpx_status status;
  size_t nodes = crpx_tree_node_count(tree);
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !nodes || !tips || !taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (nodes > SIZE_MAX/sizeof(*descendants)) return CRPX_ERROR_OVERFLOW;
  for (size_t i = 0; i < taxa; ++i)
    if (tips[i] >= nodes || crpx_tree_first_child(tree, tips[i]) != CRPX_NODE_NONE) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &w);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh(w, tree);
  if (status == CRPX_OK) status = crpx_hashtable_new(nodes, 0, &unique);
  if (status != CRPX_OK) goto done;
  s = calloc(1, sizeof(*s)); descendants = calloc(nodes, sizeof(*descendants));
  if (!s || !descendants) { status = CRPX_ERROR_ALLOCATION; goto done; }
  s->taxa = taxa;
  s->splits = calloc(nodes, sizeof(*s->splits));
  if (!s->splits) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t node = 0; node < nodes; ++node) {
    status = crpx_bipartition_new(taxa, descendants+node);
    if (status != CRPX_OK) goto done;
  }
  for (size_t i = 0; i < taxa; ++i) (void)crpx_bipartition_set(descendants[tips[i]], i, true);
  const crpx_node_id *order = crpx_tree_workspace_postorder(w);
  for (size_t i = 0; i < nodes; ++i) {
    crpx_node_id node = order[i];
    for (crpx_node_id child = crpx_tree_first_child(tree,node); child != CRPX_NODE_NONE;
         child = crpx_tree_next_sibling(tree,child))
      (void)crpx_bipartition_combine(descendants[node], descendants[node], descendants[child], CRPX_BITS_OR);
  }
  for (size_t node = 0; node < nodes; ++node) {
    size_t ones = crpx_bipartition_count(descendants[node]);
    if (ones < 2 || taxa-ones < 2) continue;
    crpx_bipartition_canonicalise(descendants[node]);
    size_t nwords;
    const uint64_t *words = crpx_bipartition_words(descendants[node], &nwords);
    uint64_t index;
    status = crpx_hashtable_get(unique, words, nwords*sizeof(uint64_t), &index);
    if (status == CRPX_OK) continue;
    status = crpx_hashtable_put(unique, words, nwords*sizeof(uint64_t), s->count, NULL);
    if (status != CRPX_OK) goto done;
    status = crpx_bipartition_clone(descendants[node], s->splits+s->count);
    if (status != CRPX_OK) goto done;
    ++s->count;
  }
done:
  if (descendants) for (size_t node = 0; node < nodes; ++node) crpx_bipartition_free(descendants[node]);
  free(descendants); crpx_hashtable_free(unique); crpx_tree_workspace_free(w);
  if (status != CRPX_OK) crpx_splitset_free(s); else *out = s;
  return status;
}

crpx_status crpx_splitset_compare(const crpx_splitset *a, const crpx_splitset *b, size_t *rf, uint64_t *assignment_cost)
{
  if (!a || !b || a->taxa != b->taxa || (!rf && !assignment_cost)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t *left = NULL, *right = NULL, *assignment = NULL;
  double *costs = NULL;
  crpx_bipartition_t difference = NULL;
  crpx_hashtable_t common = NULL;
  bool *matched = NULL;
  crpx_status status = CRPX_OK;
  if (a->count > SIZE_MAX-b->count) return CRPX_ERROR_OVERFLOW;
  if (a->count > SIZE_MAX/sizeof(size_t) || b->count > SIZE_MAX/sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  if (assignment_cost) {
    left = a->count ? malloc(a->count*sizeof(size_t)) : NULL;
    right = b->count ? malloc(b->count*sizeof(size_t)) : NULL;
    if ((a->count && !left) || (b->count && !right)) { status = CRPX_ERROR_ALLOCATION; goto done; }
  }
  matched = b->count ? calloc(b->count, sizeof(*matched)) : NULL;
  if (b->count && !matched) { status = CRPX_ERROR_ALLOCATION; goto done; }
  status = crpx_hashtable_new(b->count, 0, &common);
  if (status != CRPX_OK) goto done;
  for (size_t j = 0; j < b->count; ++j) {
    size_t nwords;
    const uint64_t *words = crpx_bipartition_words(b->splits[j], &nwords);
    status = crpx_hashtable_put(common, words, nwords*sizeof(uint64_t), j, NULL);
    if (status != CRPX_OK) goto done;
  }
  size_t nl = 0, nr = 0;
  for (size_t i = 0; i < a->count; ++i) {
    size_t nwords;
    uint64_t index;
    const uint64_t *words = crpx_bipartition_words(a->splits[i], &nwords);
    status = crpx_hashtable_get(common, words, nwords*sizeof(uint64_t), &index);
    if (status == CRPX_ERROR_NOT_FOUND) {
      if (left) left[nl] = i;
      ++nl;
    } else if (status == CRPX_OK) matched[index] = true;
    else goto done;
  }
  status = CRPX_OK;
  for (size_t j = 0; j < b->count; ++j) if (!matched[j]) {
    if (right) right[nr] = j;
    ++nr;
  }
  uint64_t hdist = 0;
  size_t n = nl > nr ? nl : nr;
  if (assignment_cost && n) {
    if (n > SIZE_MAX/n || n*n > SIZE_MAX/sizeof(double) ||
        (long double)n*(a->taxa/2) > 9007199254740992.0L) { status = CRPX_ERROR_OVERFLOW; goto done; }
    costs = calloc(n*n, sizeof(double)); assignment = malloc(n*sizeof(size_t));
    if (!costs || !assignment) { status = CRPX_ERROR_ALLOCATION; goto done; }
    status = crpx_bipartition_new(a->taxa, &difference);
    if (status != CRPX_OK) goto done;
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
      size_t cost = 0;
      if (i < nl && j < nr) {
        (void)crpx_bipartition_combine(difference, a->splits[left[i]], b->splits[right[j]], CRPX_BITS_XOR);
        cost = crpx_bipartition_count(difference);
        if (a->taxa-cost < cost) cost = a->taxa-cost;
      } else if (i < nl) cost = crpx_bipartition_count(a->splits[left[i]]);
      else if (j < nr) cost = crpx_bipartition_count(b->splits[right[j]]);
      costs[i*n+j] = (double)cost;
    }
    double total;
    status = crpx_assignment_solve(costs, n, n, assignment, &total);
    if (status != CRPX_OK) goto done;
    hdist = (uint64_t)total;
  }
  if (rf) *rf = nl+nr;
  if (assignment_cost) *assignment_cost = hdist;
done:
  free(left); free(right); free(assignment); free(costs); free(matched);
  crpx_hashtable_free(common); crpx_bipartition_free(difference);
  return status;
}

crpx_status crpx_gene_species_split_distances(const crpx_tree *gene, const crpx_tree *species,
                                             const crpx_node_id *mapping, size_t *rf, uint64_t *assignment_cost)
{
  size_t taxa = crpx_tree_leaf_count(gene), count = 0;
  if (!gene || !species || !mapping || !taxa || (!rf && !assignment_cost)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *gene_tips = malloc(taxa*sizeof(*gene_tips)), *species_tips = malloc(taxa*sizeof(*species_tips));
  crpx_splitset *a = NULL, *b = NULL;
  if (!gene_tips || !species_tips) { free(gene_tips); free(species_tips); return CRPX_ERROR_ALLOCATION; }
  for (size_t node = 0; node < crpx_tree_node_count(gene); ++node)
    if (crpx_tree_first_child(gene, (crpx_node_id)node) == CRPX_NODE_NONE) {
      gene_tips[count] = (crpx_node_id)node; species_tips[count++] = mapping[node];
    }
  crpx_status status = crpx_splitset_from_tree(gene, gene_tips, taxa, &a);
  if (status == CRPX_OK) status = crpx_splitset_from_tree(species, species_tips, taxa, &b);
  if (status == CRPX_OK) status = crpx_splitset_compare(a, b, rf, assignment_cost);
  crpx_splitset_free(a); crpx_splitset_free(b); free(gene_tips); free(species_tips);
  return status;
}
