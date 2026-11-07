/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "character_data.h"
#include "char_vector.h"
#include "hashtable.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_character_data {
  size_t taxa, sites, states, patterns, total_weight;
  crpx_char_vector_t names;
  crpx_hashtable_t name_index;
  double *values;
  size_t *weights, *site_patterns;
};
typedef enum { INPUT_STATES, INPUT_MASKS, INPUT_BINARY, INPUT_VALUES } input_type;

void crpx_character_data_free(crpx_character_data *data)
{
  if (!data) return;
  crpx_char_vector_free(data->names); crpx_hashtable_free(data->name_index);
  free(data->values); free(data->weights); free(data->site_patterns); free(data);
}

static crpx_status create_data(const char *const *names, size_t taxa, size_t sites, size_t states,
                               const void *observations, const size_t *weights, bool compact,
                               input_type type, crpx_character_data **out)
{
  crpx_character_data *data = NULL;
  crpx_hashtable_t patterns = NULL;
  double *column = NULL;
  crpx_status status = CRPX_OK;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!names || !taxa || !states || (sites && !observations) || (type == INPUT_MASKS && states > 64))
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX/states || taxa*states > SIZE_MAX/sizeof(double) ||
      (sites && taxa*states > SIZE_MAX/sites/sizeof(double)) || sites > SIZE_MAX/sizeof(size_t))
    return CRPX_ERROR_OVERFLOW;
  size_t width = taxa*states;
  data = calloc(1, sizeof(*data));
  if (!data) return CRPX_ERROR_ALLOCATION;
  data->taxa = taxa; data->sites = sites; data->states = states;
  status = crpx_char_vector_new(taxa, &data->names);
  if (status == CRPX_OK) status = crpx_hashtable_new(taxa, 0, &data->name_index);
  for (size_t taxon = 0; taxon < taxa && status == CRPX_OK; ++taxon) {
    if (!names[taxon] || !*names[taxon]) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
    for (const unsigned char *p = (const unsigned char *)names[taxon]; *p; ++p)
      if (*p < 32 || *p == 127) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
    if (status != CRPX_OK) break;
    bool replaced;
    status = crpx_hashtable_put_string(data->name_index, names[taxon], taxon, &replaced);
    if (status == CRPX_OK && replaced) status = CRPX_ERROR_INVALID_ARGUMENT;
    if (status == CRPX_OK) status = crpx_char_vector_set(data->names, taxon, names[taxon]);
  }
  if (status != CRPX_OK || !sites) goto done;
  data->values = malloc(sites*width*sizeof(double));
  data->weights = calloc(sites, sizeof(size_t));
  data->site_patterns = malloc(sites*sizeof(size_t));
  column = malloc(width*sizeof(double));
  if (!data->values || !data->weights || !data->site_patterns || !column) { status = CRPX_ERROR_ALLOCATION; goto done; }
  if (compact) status = crpx_hashtable_new(sites, 0, &patterns);
  if (status != CRPX_OK) goto done;
  for (size_t site = 0; site < sites; ++site) {
    for (size_t taxon = 0; taxon < taxa; ++taxon) {
      size_t offset = site*taxa+taxon;
      bool positive = false;
      uint64_t mask = type == INPUT_MASKS ? ((const uint64_t *)observations)[offset] : 0;
      size_t code = type == INPUT_STATES ? ((const size_t *)observations)[offset] :
                    type == INPUT_BINARY ? ((const uint8_t *)observations)[offset] : SIZE_MAX;
      if (type == INPUT_BINARY && code == UINT8_MAX) code = SIZE_MAX;
      if ((type == INPUT_STATES || type == INPUT_BINARY) && code != SIZE_MAX && code >= states) {
        status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
      }
      if (type == INPUT_MASKS && (!mask || (states < 64 && mask >> states))) {
        status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
      }
      for (size_t state = 0; state < states; ++state) {
        double value = type == INPUT_VALUES ? ((const double *)observations)[offset*states+state] :
                       type == INPUT_MASKS ? ((mask & (UINT64_C(1) << state)) ? 1 : 0) :
                       code == SIZE_MAX || code == state ? 1 : 0;
        if (!isfinite(value) || value < 0) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
        column[taxon*states+state] = value == 0 ? 0 : value;
        positive |= value > 0;
      }
      if (!positive) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    }
    uint64_t pattern = data->patterns;
    if (!compact || crpx_hashtable_get(patterns, column, width*sizeof(double), &pattern) == CRPX_ERROR_NOT_FOUND) {
      pattern = data->patterns;
      memcpy(data->values+data->patterns*width, column, width*sizeof(double));
      if (compact) status = crpx_hashtable_put(patterns, column, width*sizeof(double), pattern, NULL);
      if (status != CRPX_OK) goto done;
      ++data->patterns;
    }
    size_t weight = weights ? weights[site] : 1;
    if (data->total_weight > SIZE_MAX-weight) { status = CRPX_ERROR_OVERFLOW; goto done; }
    data->total_weight += weight;
    data->weights[pattern] += weight; /* Bounded by the checked total. */
    data->site_patterns[site] = (size_t)pattern;
  }
  {
    double *values = realloc(data->values, data->patterns*width*sizeof(double));
    size_t *pattern_weights = realloc(data->weights, data->patterns*sizeof(size_t));
    if (values) data->values = values;
    if (pattern_weights) data->weights = pattern_weights;
  }
done:
  free(column); crpx_hashtable_free(patterns);
  if (status != CRPX_OK) crpx_character_data_free(data); else *out = data;
  return status;
}

crpx_status crpx_character_data_from_states(const char *const *names, size_t taxa, size_t sites, size_t states,
                                           const size_t *values, const size_t *weights, bool compact,
                                           crpx_character_data **out)
{ return create_data(names,taxa,sites,states,values,weights,compact,INPUT_STATES,out); }
crpx_status crpx_character_data_from_masks(const char *const *names, size_t taxa, size_t sites, size_t states,
                                          const uint64_t *values, const size_t *weights, bool compact,
                                          crpx_character_data **out)
{ return create_data(names,taxa,sites,states,values,weights,compact,INPUT_MASKS,out); }
crpx_status crpx_character_data_from_binary(const char *const *names, size_t taxa, size_t sites,
                                           const uint8_t *values, const size_t *weights, bool compact,
                                           crpx_character_data **out)
{ return create_data(names,taxa,sites,2,values,weights,compact,INPUT_BINARY,out); }
crpx_status crpx_character_data_from_tip_likelihoods(const char *const *names, size_t taxa, size_t sites,
                                                   size_t states, const double *values, const size_t *weights,
                                                   bool compact, crpx_character_data **out)
{ return create_data(names,taxa,sites,states,values,weights,compact,INPUT_VALUES,out); }

crpx_status crpx_character_data_from_alignment_dna(const crpx_alignment *alignment, bool compact,
                                                 crpx_character_data **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!alignment || !crpx_alignment_is_aligned(alignment)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa = crpx_alignment_taxa(alignment), sites = crpx_alignment_sites(alignment);
  if (!taxa || taxa > SIZE_MAX/sizeof(char *) || (sites && taxa > SIZE_MAX/sites/sizeof(uint64_t)))
    return CRPX_ERROR_OVERFLOW;
  const char **names = malloc(taxa*sizeof(*names));
  uint64_t *masks = sites ? malloc(taxa*sites*sizeof(*masks)) : NULL;
  if (!names || (sites && !masks)) { free(names); free(masks); return CRPX_ERROR_ALLOCATION; }
  crpx_status status = CRPX_OK;
  for (size_t taxon = 0; taxon < taxa; ++taxon) {
    names[taxon] = crpx_alignment_name(alignment,taxon);
    const unsigned char *sequence = (const unsigned char *)crpx_alignment_sequence(alignment,taxon);
    for (size_t site = 0; site < sites; ++site) {
      masks[site*taxa+taxon] = crpx_dna_state_mask(sequence[site]);
      if (!masks[site*taxa+taxon]) status = CRPX_ERROR_INVALID_ARGUMENT;
    }
  }
  if (status == CRPX_OK) status = crpx_character_data_from_masks(names,taxa,sites,4,masks,NULL,compact,out);
  free(names); free(masks); return status;
}

size_t crpx_character_data_taxa(const crpx_character_data *data) { return data ? data->taxa : 0; }
size_t crpx_character_data_sites(const crpx_character_data *data) { return data ? data->sites : 0; }
size_t crpx_character_data_states(const crpx_character_data *data) { return data ? data->states : 0; }
size_t crpx_character_data_patterns(const crpx_character_data *data) { return data ? data->patterns : 0; }
size_t crpx_character_data_total_weight(const crpx_character_data *data) { return data ? data->total_weight : 0; }
const char *crpx_character_data_name(const crpx_character_data *data, size_t taxon)
{ return data ? crpx_char_vector_get(data->names,taxon) : NULL; }
const double *crpx_character_data_pattern_data(const crpx_character_data *data) { return data ? data->values : NULL; }
const size_t *crpx_character_data_pattern_weights(const crpx_character_data *data) { return data ? data->weights : NULL; }
const size_t *crpx_character_data_site_patterns(const crpx_character_data *data) { return data ? data->site_patterns : NULL; }
crpx_status crpx_character_data_find_taxon(const crpx_character_data *data, const char *name, size_t *taxon)
{
  if (!data || !name || !taxon) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t value;
  crpx_status status = crpx_hashtable_get_string(data->name_index,name,&value);
  if (status == CRPX_OK) *taxon = (size_t)value;
  return status;
}

crpx_status crpx_character_data_map_tips(const crpx_character_data *data, const crpx_tree *tree,
                                       crpx_node_id *taxa_to_nodes, size_t capacity)
{
  if (!data || !tree || !taxa_to_nodes) return CRPX_ERROR_INVALID_ARGUMENT;
  if (capacity < data->taxa) return CRPX_ERROR_OUT_OF_RANGE;
  if (crpx_tree_leaf_count(tree) != data->taxa) return CRPX_ERROR_TOPOLOGY;
  crpx_status status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  if (data->taxa > SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *mapping = malloc(data->taxa*sizeof(*mapping));
  if (!mapping) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < data->taxa; ++i) mapping[i] = CRPX_NODE_NONE;
  for (size_t node = 0; node < crpx_tree_node_count(tree); ++node) {
    if (crpx_tree_first_child(tree,(crpx_node_id)node) != CRPX_NODE_NONE) continue;
    size_t taxon;
    status = crpx_character_data_find_taxon(data,crpx_tree_label(tree,(crpx_node_id)node),&taxon);
    if (status != CRPX_OK || mapping[taxon] != CRPX_NODE_NONE) { status = CRPX_ERROR_TOPOLOGY; break; }
    mapping[taxon] = (crpx_node_id)node;
  }
  if (status == CRPX_OK) memcpy(taxa_to_nodes,mapping,data->taxa*sizeof(*mapping));
  free(mapping); return status;
}
