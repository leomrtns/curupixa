/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_CHARACTER_DATA_H
#define CURUPIXA_CHARACTER_DATA_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "alignment.h"
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_character_data crpx_character_data;
/* Immutable owned data with unique nonempty names and any positive state count.
 * Input observations are site-major [site*taxa+taxon]. SIZE_MAX denotes unknown
 * in from_states; binary observations accept 0,1 and UINT8_MAX (unknown).
 * Masks require states<=64 and nonzero masks; use all valid bits for unknown.
 * weights may be NULL (all one), including zero site weights is permitted.
 * Exact repeated columns optionally share a pattern and sum their checked weights.
 * No input arrays are retained. All constructor failure outputs become NULL. */
crpx_status crpx_character_data_from_states(const char *const *names, size_t taxa, size_t sites, size_t states,
                                           const size_t *observations, const size_t *weights, bool compact,
                                           crpx_character_data **out);
crpx_status crpx_character_data_from_masks(const char *const *names, size_t taxa, size_t sites, size_t states,
                                          const uint64_t *masks, const size_t *weights, bool compact,
                                          crpx_character_data **out);
crpx_status crpx_character_data_from_binary(const char *const *names, size_t taxa, size_t sites,
                                           const uint8_t *observations, const size_t *weights, bool compact,
                                           crpx_character_data **out);
/* Tip likelihoods are nonnegative finite observation weights, NOT normalised
 * probabilities: at least one state per tip/site must be positive. Input layout
 * is [site][taxon][state]. Unknown observations use all ones. Parsimony treats
 * positive entries as allowed states and ignores their magnitudes. */
crpx_status crpx_character_data_from_tip_likelihoods(const char *const *names, size_t taxa, size_t sites,
                                                   size_t states, const double *values, const size_t *weights,
                                                   bool compact, crpx_character_data **out);
crpx_status crpx_character_data_from_alignment_dna(const crpx_alignment *alignment, bool compact,
                                                 crpx_character_data **out);
void crpx_character_data_free(crpx_character_data *data);
size_t crpx_character_data_taxa(const crpx_character_data *data);
size_t crpx_character_data_sites(const crpx_character_data *data);
size_t crpx_character_data_states(const crpx_character_data *data);
size_t crpx_character_data_patterns(const crpx_character_data *data);
size_t crpx_character_data_total_weight(const crpx_character_data *data);
const char *crpx_character_data_name(const crpx_character_data *data, size_t taxon);
crpx_status crpx_character_data_find_taxon(const crpx_character_data *data, const char *name, size_t *taxon);
/* Borrowed read-only arrays, valid until free. Dense layout [pattern][taxon][state].
 * Storage is O(patterns*taxa*states), independent of tree and numerical workspaces. */
const double *crpx_character_data_pattern_data(const crpx_character_data *data);
const size_t *crpx_character_data_pattern_weights(const crpx_character_data *data);
const size_t *crpx_character_data_site_patterns(const crpx_character_data *data);
/* Map data taxa to unique, exactly matching tree leaf labels. Output remains
 * unchanged on failure. Neither tree nor output is retained by this function. */
crpx_status crpx_character_data_map_tips(const crpx_character_data *data, const crpx_tree *tree,
                                       crpx_node_id *taxa_to_nodes, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
