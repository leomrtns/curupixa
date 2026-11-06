/* Reimplementation of biomcmc-lib reconciliation (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "reconciliation.h"
#include <stdlib.h>

struct crpx_reconciliation {
  size_t gene_nodes, species_nodes;
  uint64_t duplications, losses, extra_lineages;
  crpx_node_id *mapping;
  uint8_t *duplication_nodes;
  uint64_t *node_losses, *lineages;
};

static bool is_binary(const crpx_tree *tree)
{
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    crpx_node_id a = crpx_tree_first_child(tree, (crpx_node_id)i);
    if (a == CRPX_NODE_NONE) continue;
    crpx_node_id b = crpx_tree_next_sibling(tree, a);
    if (b == CRPX_NODE_NONE || crpx_tree_next_sibling(tree, b) != CRPX_NODE_NONE) return false;
  }
  return true;
}

static crpx_node_id local_lca(const crpx_tree *species, const uint32_t *depth, crpx_node_id a, crpx_node_id b)
{
  while (depth[a] > depth[b]) a = crpx_tree_parent(species, a);
  while (depth[b] > depth[a]) b = crpx_tree_parent(species, b);
  while (a != b) { a = crpx_tree_parent(species, a); b = crpx_tree_parent(species, b); }
  return a;
}

crpx_status crpx_reconcile_lca(const crpx_tree *gene, const crpx_tree *species,
                              const crpx_node_id *leaf_map, crpx_reconciliation **out)
{
  crpx_tree_workspace *gw = NULL, *sw = NULL;
  crpx_reconciliation *r = NULL;
  crpx_status status;
  size_t gn = crpx_tree_node_count(gene), sn = crpx_tree_node_count(species);
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!gene || !species || !gn || !sn || !leaf_map) return CRPX_ERROR_INVALID_ARGUMENT;
  if (gn > SIZE_MAX/sizeof(uint64_t) || sn > SIZE_MAX/sizeof(uint64_t)) return CRPX_ERROR_OVERFLOW;
  status = crpx_tree_validate(gene);
  if (status == CRPX_OK) status = crpx_tree_validate(species);
  if (status != CRPX_OK) return status;
  if (!is_binary(gene) || !is_binary(species)) return CRPX_ERROR_UNSUPPORTED;
  for (size_t node = 0; node < gn; ++node) if (crpx_tree_first_child(gene, (crpx_node_id)node) == CRPX_NODE_NONE) {
    crpx_node_id mapped = leaf_map[node];
    if (mapped >= sn || crpx_tree_first_child(species, mapped) != CRPX_NODE_NONE) return CRPX_ERROR_INVALID_ARGUMENT;
  }
  r = calloc(1, sizeof(*r));
  if (!r) return CRPX_ERROR_ALLOCATION;
  r->gene_nodes = gn; r->species_nodes = sn;
  r->mapping = malloc(gn*sizeof(*r->mapping));
  r->duplication_nodes = calloc(gn, sizeof(*r->duplication_nodes));
  r->node_losses = calloc(gn, sizeof(*r->node_losses));
  r->lineages = calloc(sn, sizeof(*r->lineages));
  if (!r->mapping || !r->duplication_nodes || !r->node_losses || !r->lineages) { status = CRPX_ERROR_ALLOCATION; goto done; }
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &gw);
  if (status == CRPX_OK) status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_DEPTH, &sw);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh(gw, gene);
  if (status == CRPX_OK) status = crpx_tree_workspace_refresh(sw, species);
  if (status != CRPX_OK) goto done;
  const crpx_node_id *postorder = crpx_tree_workspace_postorder(gw);
  const uint32_t *depth = crpx_tree_workspace_depth(sw);
  for (size_t i = 0; i < gn; ++i) {
    crpx_node_id node = postorder[i], left = crpx_tree_first_child(gene, node);
    if (left == CRPX_NODE_NONE) { r->mapping[node] = leaf_map[node]; continue; }
    crpx_node_id right = crpx_tree_next_sibling(gene, left);
    crpx_node_id a = r->mapping[left], b = r->mapping[right];
    crpx_node_id parent = local_lca(species, depth, a, b);
    r->mapping[node] = parent;
    bool duplication = parent == a || parent == b;
    r->duplication_nodes[node] = duplication;
    r->duplications += duplication;
    uint64_t paths = (uint64_t)(depth[a]-depth[parent]) + (uint64_t)(depth[b]-depth[parent]);
    uint64_t losses = paths - (duplication ? 0 : 2);
    if (r->losses > UINT64_MAX-losses) { status = CRPX_ERROR_OVERFLOW; goto done; }
    r->node_losses[node] = losses;
    r->losses += losses;
  }
  for (size_t node = 0; node < gn; ++node) {
    crpx_node_id parent = crpx_tree_parent(gene, (crpx_node_id)node);
    if (parent == CRPX_NODE_NONE) continue;
    crpx_node_id to = r->mapping[parent];
    for (crpx_node_id from = r->mapping[node]; from != to; from = crpx_tree_parent(species, from)) ++r->lineages[from];
  }
  for (size_t node = 0; node < sn; ++node) if (r->lineages[node] > 1) {
    uint64_t extra = r->lineages[node]-1;
    if (r->extra_lineages > UINT64_MAX-extra) { status = CRPX_ERROR_OVERFLOW; goto done; }
    r->extra_lineages += extra;
  }
done:
  crpx_tree_workspace_free(gw); crpx_tree_workspace_free(sw);
  if (status != CRPX_OK) crpx_reconciliation_free(r);
  else *out = r;
  return status;
}

void crpx_reconciliation_free(crpx_reconciliation *r)
{
  if (r) { free(r->mapping); free(r->duplication_nodes); free(r->node_losses); free(r->lineages); free(r); }
}
uint64_t crpx_reconciliation_duplications(const crpx_reconciliation *r) { return r ? r->duplications : 0; }
uint64_t crpx_reconciliation_losses(const crpx_reconciliation *r) { return r ? r->losses : 0; }
uint64_t crpx_reconciliation_extra_lineages(const crpx_reconciliation *r) { return r ? r->extra_lineages : 0; }
size_t crpx_reconciliation_gene_nodes(const crpx_reconciliation *r) { return r ? r->gene_nodes : 0; }
const crpx_node_id *crpx_reconciliation_mapping(const crpx_reconciliation *r) { return r ? r->mapping : NULL; }
const uint8_t *crpx_reconciliation_duplication_nodes(const crpx_reconciliation *r) { return r ? r->duplication_nodes : NULL; }
const uint64_t *crpx_reconciliation_node_losses(const crpx_reconciliation *r) { return r ? r->node_losses : NULL; }
const uint64_t *crpx_reconciliation_species_lineages(const crpx_reconciliation *r, size_t *nodes)
{
  if (nodes) *nodes = r ? r->species_nodes : 0;
  return r ? r->lineages : NULL;
}
