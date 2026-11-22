/* Two-direction reconciliation derived from guenomu/genefam-dist topology_mrca.
 * Copyright (C) 2009-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Score alternative gene roots and extract their mappings using reusable directional messages. */
#include "reconciliation_roots.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
  crpx_node_id species_node;
  uint64_t duplications, losses, represented_losses;
} reconciliation_message;

typedef struct {
  crpx_node_id neighbours[3];
  size_t degree;
  reconciliation_message outgoing[3];
} reconciliation_vertex;

struct crpx_reconciliation_workspace {
  reconciliation_vertex *vertices;
  crpx_node_id *gene_parent, *gene_order, *species_order;
  uint64_t *presence;
  uint32_t *represented_depth;
  size_t *depth_counts;
  uint8_t *represented_children;
  size_t gene_capacity, species_capacity, gene_nodes;
  crpx_node_id artificial_root;
  crpx_root_reconciliations result;
  bool valid;
};

static bool
add_counts (uint64_t first, uint64_t second, uint64_t *sum)
{
  if (first > UINT64_MAX - second) return false;
  *sum = first + second;
  return true;
}

/* Combine two sides of a gene branching. Speciations save the two outgoing species edges; duplications retain their
   full path lengths as losses. */
static crpx_status
combine_sides (const crpx_ancestry_index *species, const uint32_t *represented_depth,
               reconciliation_message first, reconciliation_message second, reconciliation_message *result)
{
  crpx_node_id ancestor;
  crpx_status status = crpx_ancestry_index_lca (species, first.species_node, second.species_node, &ancestor);
  if (status != CRPX_OK) return status;
  bool duplication = ancestor == first.species_node || ancestor == second.species_node;
  uint64_t full_path = (uint64_t) (crpx_ancestry_index_depth (species, first.species_node) -
                                  crpx_ancestry_index_depth (species, ancestor)) +
    crpx_ancestry_index_depth (species, second.species_node) - crpx_ancestry_index_depth (species, ancestor);
  uint64_t represented_path = (uint64_t) (represented_depth[first.species_node] - represented_depth[ancestor]) +
    represented_depth[second.species_node] - represented_depth[ancestor];
  if (!duplication && (full_path < 2 || represented_path < 2)) return CRPX_ERROR_TOPOLOGY;
  reconciliation_message combined = {ancestor, 0, 0, 0};
  uint64_t speciation_edges = duplication ? 0 : 2;
  if (!add_counts (first.duplications, second.duplications, &combined.duplications) ||
      !add_counts (combined.duplications, duplication, &combined.duplications) ||
      !add_counts (first.losses, second.losses, &combined.losses) ||
      !add_counts (combined.losses, full_path - speciation_edges, &combined.losses) ||
      !add_counts (first.represented_losses, second.represented_losses, &combined.represented_losses) ||
      !add_counts (combined.represented_losses, represented_path - speciation_edges, &combined.represented_losses))
    return CRPX_ERROR_OVERFLOW;
  *result = combined;
  return CRPX_OK;
}

static size_t
neighbour_slot (const reconciliation_vertex *vertex, crpx_node_id neighbour)
{
  for (size_t slot = 0; slot < vertex->degree; slot++)
    if (vertex->neighbours[slot] == neighbour) return slot;
  return SIZE_MAX;
}

static crpx_status
connect_vertices (reconciliation_vertex *vertices, crpx_node_id first, crpx_node_id second)
{
  if (vertices[first].degree == 3 || vertices[second].degree == 3) return CRPX_ERROR_UNSUPPORTED;
  vertices[first].neighbours[vertices[first].degree++] = second;
  vertices[second].neighbours[vertices[second].degree++] = first;
  return CRPX_OK;
}

static crpx_status
combine_component (reconciliation_vertex *vertices, crpx_node_id node, size_t excluded_slot,
                   const crpx_node_id *tip_map, const crpx_ancestry_index *species,
                   const uint32_t *represented_depth, reconciliation_message *result)
{
  reconciliation_vertex *vertex = vertices + node;
  if (vertex->degree == 1) {
    *result = (reconciliation_message) {tip_map[node], 0, 0, 0};
    return CRPX_OK;
  }
  if (vertex->degree != 3) return CRPX_ERROR_TOPOLOGY;
  reconciliation_message sides[2];
  size_t side_count = 0;
  for (size_t slot = 0; slot < vertex->degree; slot++) if (slot != excluded_slot) {
    crpx_node_id neighbour = vertex->neighbours[slot];
    size_t reverse_slot = neighbour_slot (vertices + neighbour, node);
    sides[side_count++] = vertices[neighbour].outgoing[reverse_slot];
  }
  return combine_sides (species, represented_depth, sides[0], sides[1], result);
}

crpx_status
crpx_reconciliation_workspace_new (crpx_reconciliation_workspace **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = calloc (1, sizeof (**result));
  return *result ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}

void
crpx_reconciliation_workspace_clear (crpx_reconciliation_workspace *workspace)
{
  if (!workspace) return;
  free (workspace->vertices);
  free (workspace->gene_parent);
  free (workspace->gene_order);
  free (workspace->species_order);
  free (workspace->presence);
  free (workspace->represented_depth);
  free (workspace->depth_counts);
  free (workspace->represented_children);
  free (workspace->result.scores);
  memset (workspace, 0, sizeof (*workspace));
}

void
crpx_reconciliation_workspace_free (crpx_reconciliation_workspace *workspace)
{
  if (!workspace) return;
  crpx_reconciliation_workspace_clear (workspace);
  free (workspace);
}

size_t
crpx_reconciliation_workspace_storage_size (const crpx_reconciliation_workspace *workspace)
{
  if (!workspace) return 0;
  return sizeof (*workspace) + workspace->gene_capacity *
    (sizeof (reconciliation_vertex) + 2 * sizeof (crpx_node_id) + sizeof (crpx_root_reconciliation_score)) +
    workspace->species_capacity * (sizeof (crpx_node_id) + sizeof (uint64_t) + sizeof (uint32_t) +
                                  sizeof (size_t) + sizeof (uint8_t)) +
    (workspace->species_capacity ? sizeof (size_t) : 0);
}

static crpx_status
reserve_workspace (crpx_reconciliation_workspace *workspace, size_t gene_nodes, size_t species_nodes)
{
  if (gene_nodes > SIZE_MAX / sizeof (reconciliation_vertex) ||
      gene_nodes > SIZE_MAX / sizeof (crpx_root_reconciliation_score) ||
      species_nodes == SIZE_MAX || species_nodes + 1 > SIZE_MAX / sizeof (size_t) ||
      species_nodes > SIZE_MAX / sizeof (uint64_t)) return CRPX_ERROR_OVERFLOW;
  if (gene_nodes > workspace->gene_capacity) {
    reconciliation_vertex *vertices = malloc (gene_nodes * sizeof (*vertices));
    crpx_node_id *parent = malloc (gene_nodes * sizeof (*parent));
    crpx_node_id *order = malloc (gene_nodes * sizeof (*order));
    crpx_root_reconciliation_score *scores = malloc (gene_nodes * sizeof (*scores));
    if (!vertices || !parent || !order || !scores) {
      free (vertices); free (parent); free (order); free (scores);
      return CRPX_ERROR_ALLOCATION;
    }
    free (workspace->vertices); free (workspace->gene_parent); free (workspace->gene_order);
    free (workspace->result.scores);
    workspace->vertices = vertices;
    workspace->gene_parent = parent;
    workspace->gene_order = order;
    workspace->result.scores = scores;
    workspace->gene_capacity = gene_nodes;
  }
  if (species_nodes > workspace->species_capacity) {
    crpx_node_id *order = malloc (species_nodes * sizeof (*order));
    uint64_t *presence = malloc (species_nodes * sizeof (*presence));
    uint32_t *depth = malloc (species_nodes * sizeof (*depth));
    size_t *counts = malloc ((species_nodes + 1) * sizeof (*counts));
    uint8_t *children = malloc (species_nodes * sizeof (*children));
    if (!order || !presence || !depth || !counts || !children) {
      free (order); free (presence); free (depth); free (counts); free (children);
      return CRPX_ERROR_ALLOCATION;
    }
    free (workspace->species_order); free (workspace->presence); free (workspace->represented_depth);
    free (workspace->depth_counts); free (workspace->represented_children);
    workspace->species_order = order;
    workspace->presence = presence;
    workspace->represented_depth = depth;
    workspace->depth_counts = counts;
    workspace->represented_children = children;
    workspace->species_capacity = species_nodes;
  }
  memset (workspace->vertices, 0, gene_nodes * sizeof (*workspace->vertices));
  memset (workspace->presence, 0, species_nodes * sizeof (*workspace->presence));
  memset (workspace->depth_counts, 0, (species_nodes + 1) * sizeof (*workspace->depth_counts));
  memset (workspace->represented_children, 0, species_nodes * sizeof (*workspace->represented_children));
  return CRPX_OK;
}

/* Species nodes are sorted from shallow to deep once per evaluation. Restricting depths to represented branches
   reproduces the legacy induced-tree losses. */
static void
prepare_species_depths (crpx_reconciliation_workspace *workspace, const crpx_ancestry_index *species)
{
  size_t species_nodes = crpx_ancestry_index_nodes (species);
  size_t *counts = workspace->depth_counts;
  crpx_node_id *order = workspace->species_order;
  for (size_t node = 0; node < species_nodes; node++)
    counts[crpx_ancestry_index_depth (species, (crpx_node_id) node) + 1]++;
  for (size_t depth = 1; depth <= species_nodes; depth++) counts[depth] += counts[depth - 1];
  for (size_t node = 0; node < species_nodes; node++)
    order[counts[crpx_ancestry_index_depth (species, (crpx_node_id) node)]++] = (crpx_node_id) node;
  for (size_t position = species_nodes; position > 0;) {
    crpx_node_id node = order[--position], parent = crpx_ancestry_index_parent (species, node);
    if (parent != CRPX_NODE_NONE) workspace->presence[parent] += workspace->presence[node];
  }
  for (size_t node = 0; node < species_nodes; node++) if (workspace->presence[node]) {
    crpx_node_id parent = crpx_ancestry_index_parent (species, (crpx_node_id) node);
    if (parent != CRPX_NODE_NONE) workspace->represented_children[parent]++;
  }
  for (size_t position = 0; position < species_nodes; position++) {
    crpx_node_id node = order[position], parent = crpx_ancestry_index_parent (species, node);
    workspace->represented_depth[node] = (parent == CRPX_NODE_NONE ? 0 : workspace->represented_depth[parent]) +
      (workspace->presence[node] &&
       (crpx_ancestry_index_is_leaf (species, node) || workspace->represented_children[node] == 2));
  }
}

static crpx_status
prepare_gene_graph (crpx_reconciliation_workspace *workspace, const crpx_tree *gene)
{
  reconciliation_vertex *vertices = workspace->vertices;
  crpx_node_id root = crpx_tree_root (gene), first = crpx_tree_first_child (gene, root);
  crpx_node_id second = crpx_tree_next_sibling (gene, first);
  for (size_t node = 0; node < workspace->gene_nodes; node++) if (node != root) {
    crpx_node_id parent = crpx_tree_parent (gene, (crpx_node_id) node);
    if (parent != root) {
      crpx_status status = connect_vertices (vertices, (crpx_node_id) node, parent);
      if (status != CRPX_OK) return status;
    }
  }
  crpx_status status = connect_vertices (vertices, first, second);
  if (status != CRPX_OK) return status;
  /* Sorting at most three neighbours makes candidate order independent of child order. */
  for (size_t node = 0; node < workspace->gene_nodes; node++) {
    reconciliation_vertex *vertex = vertices + node;
    for (size_t slot = 1; slot < vertex->degree; slot++) {
      crpx_node_id neighbour = vertex->neighbours[slot];
      size_t position = slot;
      while (position && vertex->neighbours[position - 1] > neighbour) {
        vertex->neighbours[position] = vertex->neighbours[position - 1];
        position--;
      }
      vertex->neighbours[position] = neighbour;
    }
    workspace->gene_parent[node] = CRPX_NODE_NONE;
  }
  workspace->gene_order[0] = first;
  size_t visited = 1;
  for (size_t position = 0; position < visited; position++) {
    crpx_node_id node = workspace->gene_order[position];
    reconciliation_vertex *vertex = vertices + node;
    for (size_t slot = 0; slot < vertex->degree; slot++) {
      crpx_node_id neighbour = vertex->neighbours[slot];
      if (neighbour == workspace->gene_parent[node]) continue;
      workspace->gene_parent[neighbour] = node;
      workspace->gene_order[visited++] = neighbour;
    }
  }
  return visited == workspace->gene_nodes - 1 ? CRPX_OK : CRPX_ERROR_TOPOLOGY;
}

static crpx_status
calculate_directional_messages (crpx_reconciliation_workspace *workspace, const crpx_ancestry_index *species,
                                const crpx_node_id *tip_map)
{
  reconciliation_vertex *vertices = workspace->vertices;
  size_t visited = workspace->gene_nodes - 1;
  /* First accumulate towards an arbitrary temporary root, then send the complementary messages away from it. Each
   directed edge is evaluated once. */
  for (size_t position = visited; position > 1;) {
    crpx_node_id node = workspace->gene_order[--position];
    size_t slot = neighbour_slot (vertices + node, workspace->gene_parent[node]);
    crpx_status status = combine_component (vertices, node, slot, tip_map, species,
                                           workspace->represented_depth, &vertices[node].outgoing[slot]);
    if (status != CRPX_OK) return status;
  }
  for (size_t position = 0; position < visited; position++) {
    crpx_node_id node = workspace->gene_order[position];
    reconciliation_vertex *vertex = vertices + node;
    for (size_t slot = 0; slot < vertex->degree; slot++)
      if (vertex->neighbours[slot] != workspace->gene_parent[node]) {
        crpx_status status = combine_component (vertices, node, slot, tip_map, species,
                                               workspace->represented_depth, &vertex->outgoing[slot]);
        if (status != CRPX_OK) return status;
      }
  }
  return CRPX_OK;
}

crpx_status
crpx_reconcile_all_roots_into (const crpx_tree *gene, const crpx_ancestry_index *species,
                              const crpx_node_id *tip_map, crpx_reconciliation_scope scope,
                              crpx_reconciliation_workspace *workspace, const crpx_root_reconciliations **result)
{ /* Reuse directional messages to score every gene-root edge without rerooting the tree. */
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (workspace) { workspace->valid = false; workspace->result.count = 0; }
  if (!gene || !species || !tip_map || !workspace ||
      (scope != CRPX_RECONCILE_ALL_SPECIES && scope != CRPX_RECONCILE_REPRESENTED_SPECIES))
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t gene_nodes = crpx_tree_node_count (gene), species_nodes = crpx_ancestry_index_nodes (species);
  if (!gene_nodes || !crpx_ancestry_index_is_binary (species)) return CRPX_ERROR_UNSUPPORTED;
  crpx_status status = crpx_tree_validate (gene);
  if (status != CRPX_OK) return status;
  status = reserve_workspace (workspace, gene_nodes, species_nodes);
  if (status != CRPX_OK) return status;
  workspace->gene_nodes = gene_nodes;
  workspace->artificial_root = crpx_tree_root (gene);
  size_t gene_leaves = 0, represented_species = 0;
  for (size_t node = 0; node < gene_nodes; node++) {
    crpx_node_id first = crpx_tree_first_child (gene, (crpx_node_id) node);
    if (first == CRPX_NODE_NONE) {
      if (!crpx_ancestry_index_is_leaf (species, tip_map[node])) return CRPX_ERROR_INVALID_ARGUMENT;
      if (!workspace->presence[tip_map[node]]) represented_species++;
      workspace->presence[tip_map[node]]++;
      gene_leaves++;
    } else {
      crpx_node_id second = crpx_tree_next_sibling (gene, first);
      if (second == CRPX_NODE_NONE || crpx_tree_next_sibling (gene, second) != CRPX_NODE_NONE)
        return CRPX_ERROR_UNSUPPORTED;
    }
  }
  if (gene_leaves == 1) {
    workspace->vertices[workspace->artificial_root].outgoing[0].species_node = tip_map[workspace->artificial_root];
    workspace->result.scores[0] = (crpx_root_reconciliation_score) {CRPX_NODE_NONE, CRPX_NODE_NONE, 0, 0, 0};
    workspace->result.count = 1;
  } else {
    prepare_species_depths (workspace, species);
    status = prepare_gene_graph (workspace, gene);
    if (status == CRPX_OK) status = calculate_directional_messages (workspace, species, tip_map);
    if (status != CRPX_OK) return status;
    for (size_t node = 0; node < gene_nodes; node++) {
      reconciliation_vertex *vertex = workspace->vertices + node;
      for (size_t slot = 0; slot < vertex->degree; slot++) {
        crpx_node_id neighbour = vertex->neighbours[slot];
        if (node >= neighbour) continue;
        reconciliation_message combined;
        size_t reverse = neighbour_slot (workspace->vertices + neighbour, (crpx_node_id) node);
        status = combine_sides (species, workspace->represented_depth, vertex->outgoing[slot],
                                workspace->vertices[neighbour].outgoing[reverse], &combined);
        if (status != CRPX_OK) return status;
        uint64_t corrected_losses;
        if (!add_counts (combined.represented_losses, 2 * (uint64_t) (gene_leaves - represented_species),
                         &corrected_losses) || combined.duplications > UINT64_MAX / 2 ||
            corrected_losses < 2 * combined.duplications) return CRPX_ERROR_OVERFLOW;
        workspace->vertices[workspace->artificial_root].outgoing[0].species_node = combined.species_node;
        workspace->result.scores[workspace->result.count++] = (crpx_root_reconciliation_score) {
          (crpx_node_id) node, neighbour, combined.duplications,
          scope == CRPX_RECONCILE_ALL_SPECIES ? combined.losses : combined.represented_losses,
          corrected_losses - 2 * combined.duplications};
      }
    }
  }
  workspace->valid = true;
  *result = &workspace->result;
  return CRPX_OK;
}

crpx_status
crpx_reconcile_all_roots (const crpx_tree *gene, const crpx_ancestry_index *species,
                          const crpx_node_id *tip_map, crpx_reconciliation_scope scope,
                          crpx_root_reconciliations **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  crpx_reconciliation_workspace *workspace = NULL;
  const crpx_root_reconciliations *view;
  crpx_status status = crpx_reconciliation_workspace_new (&workspace);
  if (status == CRPX_OK) status = crpx_reconcile_all_roots_into (gene, species, tip_map, scope, workspace, &view);
  if (status == CRPX_OK) {
    crpx_root_reconciliations *owned = malloc (sizeof (*owned));
    if (!owned) status = CRPX_ERROR_ALLOCATION;
    else {
      *owned = workspace->result;
      workspace->result.scores = NULL;
      *result = owned;
    }
  }
  crpx_reconciliation_workspace_free (workspace);
  return status;
}

void
crpx_root_reconciliations_free (crpx_root_reconciliations *result)
{
  if (!result) return;
  free (result->scores);
  free (result);
}

const crpx_root_reconciliation_score *
crpx_root_reconciliations_scores (const crpx_root_reconciliations *result, size_t *count)
{
  if (count) *count = result ? result->count : 0;
  return result ? result->scores : NULL;
}

crpx_status
crpx_root_reconciliations_minima (const crpx_root_reconciliations *result, size_t *duplications,
                                 size_t *losses, size_t *extra_lineages)
{
  if (!result || !duplications || !losses || !extra_lineages || !result->count)
    return CRPX_ERROR_INVALID_ARGUMENT;
  size_t duplication_minimum = 0, loss_minimum = 0, lineage_minimum = 0;
  for (size_t index = 1; index < result->count; index++) {
    if (result->scores[index].duplications < result->scores[duplication_minimum].duplications)
      duplication_minimum = index;
    if (result->scores[index].losses < result->scores[loss_minimum].losses) loss_minimum = index;
    if (result->scores[index].extra_lineages < result->scores[lineage_minimum].extra_lineages)
      lineage_minimum = index;
  }
  *duplications = duplication_minimum;
  *losses = loss_minimum;
  *extra_lineages = lineage_minimum;
  return CRPX_OK;
}

crpx_status
crpx_reconciliation_root_mapping (crpx_reconciliation_workspace *workspace, crpx_node_id first,
                                  crpx_node_id second, crpx_node_id *mapping, size_t capacity)
{ /* Orient the stored messages from the selected root edge to recover all node mappings. */
  if (!workspace || !workspace->valid || !mapping) return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity < workspace->gene_nodes) return CRPX_ERROR_OUT_OF_RANGE;
  if (workspace->gene_nodes == 1) {
    if (first != CRPX_NODE_NONE || second != CRPX_NODE_NONE) return CRPX_ERROR_INVALID_ARGUMENT;
    mapping[workspace->artificial_root] = workspace->vertices[workspace->artificial_root].outgoing[0].species_node;
    return CRPX_OK;
  }
  if (first >= workspace->gene_nodes || second >= workspace->gene_nodes) return CRPX_ERROR_OUT_OF_RANGE;
  size_t first_slot = neighbour_slot (workspace->vertices + first, second);
  size_t second_slot = neighbour_slot (workspace->vertices + second, first);
  if (first_slot == SIZE_MAX || second_slot == SIZE_MAX) return CRPX_ERROR_INVALID_ARGUMENT;
  /* The artificial root is isolated from the unrooted graph; its spare message records the species MRCA, which is
   independent of the chosen gene root. */
  mapping[workspace->artificial_root] = workspace->vertices[workspace->artificial_root].outgoing[0].species_node;
  size_t visited = 2;
  workspace->gene_order[0] = first;
  workspace->gene_order[1] = second;
  workspace->gene_parent[first] = second;
  workspace->gene_parent[second] = first;
  for (size_t position = 0; position < visited; position++) {
    crpx_node_id node = workspace->gene_order[position];
    reconciliation_vertex *vertex = workspace->vertices + node;
    size_t excluded = neighbour_slot (vertex, workspace->gene_parent[node]);
    mapping[node] = vertex->outgoing[excluded].species_node;
    for (size_t slot = 0; slot < vertex->degree; slot++) {
      crpx_node_id neighbour = vertex->neighbours[slot];
      if (neighbour == workspace->gene_parent[node]) continue;
      workspace->gene_parent[neighbour] = node;
      workspace->gene_order[visited++] = neighbour;
    }
  }
  return CRPX_OK;
}
