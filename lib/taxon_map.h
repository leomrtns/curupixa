/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_TAXON_MAP_H
#define CURUPIXA_TAXON_MAP_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Exact, case-sensitive species labels, optionally extracting zero-based field
 * from a gene label separated by one delimiter. delimiter=0 requires field=0
 * and matches the complete name. Empty/missing fields, duplicate species names
 * and unmatched genes return an error. Internal map entries become NODE_NONE.
 * No substring guesses. Transactional output of gene_node_count entries. */
crpx_status crpx_gene_species_map_labels(const crpx_tree *gene,const crpx_tree *species,
                                        char delimiter,size_t field,crpx_node_id *map,size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
