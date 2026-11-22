/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: MIT */
/* Zhang et al. (2020), ASTRAL-Pro, doi:10.1093/molbev/msaa139, Definitions 1–5. Rooting and event tags are supplied by
   the application. */
/* Prepare and score ASTRAL-Pro quartet classes from rooted, tagged multi-copy gene trees. */
#include "paralog_quartets.h"
#include "ancestry_index.h"
#include <stdlib.h>

void
crpx_paralog_quartets_free (crpx_paralog_quartets *quartets)
{
  if (!quartets) return;
  free (quartets->classes);
  free (quartets);
}

static int
compare_classes (const void *left, const void *right)
{
  const crpx_paralog_quartet_class *first = left, *second = right;
  for (size_t i = 0; i < 4; i++) if (first->species[i] != second->species[i])
    return first->species[i] < second->species[i] ? -1 : 1;
  return first->anchor < second->anchor ? -1 : first->anchor != second->anchor;
}

crpx_status
crpx_paralog_quartets_new (const crpx_tree *gene, const size_t *species_of_node,
                           size_t species_taxa, const uint8_t *duplication_nodes,
                           crpx_paralog_quartets **out)
{ /* Collect speciation-driven quartets and merge equivalent gene-copy records (Zhang et al. 2020). */
  crpx_ancestry_index *ancestry = NULL;
  crpx_quartet_index *paths = NULL;
  crpx_paralog_quartets *quartets = NULL;
  crpx_node_id *tips = NULL;
  crpx_status status;
  size_t nodes = crpx_tree_node_count (gene), tip_count = crpx_tree_leaf_count (gene), capacity = 0;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!gene || !species_of_node || !duplication_nodes || !species_taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (tip_count > SIZE_MAX / sizeof (*tips)) return CRPX_ERROR_OVERFLOW;
  status = crpx_ancestry_index_new (gene, &ancestry);
  if (status != CRPX_OK) return status;
  if (!crpx_ancestry_index_is_binary (ancestry)) { status = CRPX_ERROR_UNSUPPORTED; goto done; }
  tips = malloc (tip_count * sizeof (*tips));
  quartets = calloc (1, sizeof (*quartets));
  if (!tips || !quartets) { status = CRPX_ERROR_ALLOCATION; goto done; }
  quartets->species_taxa = species_taxa;
  size_t used = 0;
  for (size_t node = 0; node < nodes; node++) {
    if (crpx_ancestry_index_is_leaf (ancestry, (crpx_node_id) node)) {
      if (species_of_node[node] >= species_taxa) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      tips[used++] = (crpx_node_id) node;
    } else if (duplication_nodes[node] > 1) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
  }
  /* A shared species below both children makes their LCA an invalid speciation. */
  for (size_t i = 0; i < tip_count; i++) for (size_t j = 0; j < i; j++)
    if (species_of_node[tips[i]] == species_of_node[tips[j]]) {
      crpx_node_id ancestor;
      status = crpx_ancestry_index_lca (ancestry, tips[i], tips[j], &ancestor);
      if (status != CRPX_OK) goto done;
      if (!duplication_nodes[ancestor]) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    }
  status = crpx_quartet_index_new (gene, tips, tip_count, &paths);
  if (status != CRPX_OK) goto done;
  for (size_t a = 0; a < tip_count; a++) for (size_t b = a + 1; b < tip_count; b++)
    for (size_t c = b + 1; c < tip_count; c++) for (size_t d = c + 1; d < tip_count; d++) {
      size_t positions[4] = {a, b, c, d};
      /* Canonical species order also fixes the topology code's interpretation. */
      for (size_t i = 1; i < 4; i++) {
        size_t position = positions[i], j = i;
        while (j && species_of_node[tips[positions[j - 1]]] > species_of_node[tips[position]]) {
          positions[j] = positions[j - 1]; j--;
        }
        positions[j] = position;
      }
      crpx_paralog_quartet_class record = {{0}, CRPX_NODE_NONE, 0};
      bool distinct = true;
      for (size_t i = 0; i < 4; i++) {
        record.species[i] = species_of_node[tips[positions[i]]];
        if (i && record.species[i] == record.species[i - 1]) distinct = false;
      }
      if (!distinct) continue;
      bool eligible = true;
      uint32_t deepest = 0;
      /* In the balanced case all triple LCAs coincide. In the unbalanced case the deepest triple excludes the outgroup and
   identifies anchor LCA. */
      for (size_t omitted = 0; omitted < 4; omitted++) {
        crpx_node_id ancestor = CRPX_NODE_NONE;
        for (size_t i = 0; i < 4; i++) if (i != omitted) {
          if (ancestor == CRPX_NODE_NONE) ancestor = tips[positions[i]];
          else {
            status = crpx_ancestry_index_lca (ancestry, ancestor, tips[positions[i]], &ancestor);
            if (status != CRPX_OK) goto done;
          }
        }
        if (duplication_nodes[ancestor]) { eligible = false; break; }
        uint32_t depth = crpx_ancestry_index_depth (ancestry, ancestor);
        if (record.anchor == CRPX_NODE_NONE || depth > deepest) { record.anchor = ancestor; deepest = depth; }
      }
      if (!eligible) continue;
      status = crpx_quartet_topology (paths, positions[0], positions[1], positions[2], positions[3], &record.topology);
      if (status != CRPX_OK) goto done;
      if (quartets->count == capacity) {
        size_t next = capacity ? capacity * 2 : 16;
        if (next < capacity || next > SIZE_MAX / sizeof (*quartets->classes)) {
          status = CRPX_ERROR_OVERFLOW; goto done;
        }
        crpx_paralog_quartet_class *grown = realloc (quartets->classes, next * sizeof (*grown));
        if (!grown) { status = CRPX_ERROR_ALLOCATION; goto done; }
        quartets->classes = grown; capacity = next;
      }
      quartets->classes[quartets->count++] = record;
    }
  if (quartets->count > 1) qsort (quartets->classes, quartets->count, sizeof (*quartets->classes), compare_classes);
  used = 0;
  for (size_t i = 0; i < quartets->count; i++) {
    if (used && !compare_classes (quartets->classes + used - 1, quartets->classes + i)) {
      if (quartets->classes[used - 1].topology != quartets->classes[i].topology) {
        status = CRPX_ERROR_TOPOLOGY; goto done;
      }
    } else quartets->classes[used++] = quartets->classes[i];
  }
  quartets->count = used;
  *out = quartets; quartets = NULL;
done:
  free (tips);
  crpx_quartet_index_free (paths);
  crpx_ancestry_index_free (ancestry);
  crpx_paralog_quartets_free (quartets);
  return status;
}

crpx_status
crpx_paralog_quartets_score (const crpx_paralog_quartets *quartets, const crpx_quartet_index *species,
                             crpx_paralog_quartet_score *score)
{ /* Score a candidate species topology against the classes prepared for this gene family. */
  crpx_paralog_quartet_score result = {0};
  if (!quartets || !species || !score || quartets->species_taxa != species->taxa)
    return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < quartets->count; i++) {
    const crpx_paralog_quartet_class *record = quartets->classes + i;
    unsigned topology;
    crpx_status status = crpx_quartet_topology (species, record->species[0], record->species[1],
                                               record->species[2], record->species[3], &topology);
    if (status != CRPX_OK) return status;
    if (result.total == UINT64_MAX) return CRPX_ERROR_OVERFLOW;
    result.total++;
    if (!topology) result.unresolved_species++;
    else if (topology == record->topology) result.agreement++;
    else result.conflict++;
  }
  *score = result;
  return CRPX_OK;
}
