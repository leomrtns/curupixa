/* Reimplementation of biomcmc-lib alignment (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_ALIGNMENT_H
#define CURUPIXA_ALIGNMENT_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "distance_matrix.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_alignment crpx_alignment;
typedef enum { CRPX_DISTANCE_P, CRPX_DISTANCE_JC69, CRPX_DISTANCE_K2P } crpx_sequence_distance;

/* Immutable owned copies, unique nonempty taxon names; ASCII sequence symbols
 * are uppercased. Ragged FASTA is allowed unless compact_patterns=true. */
crpx_status crpx_alignment_new(const char *const *names, const char *const *sequences, size_t taxa,
                               bool compact_patterns, crpx_alignment **out);
void crpx_alignment_free(crpx_alignment *alignment);
/* FASTA (including leading comment lines) or one NEXUS DATA/CHARACTERS MATRIX.
 * NEXUS accepts DNA/RNA/NUCLEOTIDE data types, quoted labels, nested comments,
 * sequential and labelled interleaved rows, standard DNA gap/missing/matchchar.
 * Blocks must be closed. Unsupported layout/format commands and multiple matrices fail.
 * Other NEXUS commands/formats return UNSUPPORTED where detected. */
crpx_status crpx_alignment_from_text(const char *text, bool compact_patterns, crpx_alignment **out);
crpx_status crpx_alignment_read(const char *path, bool compact_patterns, crpx_alignment **out);
crpx_status crpx_alignment_write_fasta(const crpx_alignment *alignment, const char *path, size_t line_width);
size_t crpx_alignment_taxa(const crpx_alignment *alignment);
bool crpx_alignment_is_aligned(const crpx_alignment *alignment);
size_t crpx_alignment_sites(const crpx_alignment *alignment);
const char *crpx_alignment_name(const crpx_alignment *alignment, size_t taxon);
const char *crpx_alignment_sequence(const crpx_alignment *alignment, size_t taxon);
crpx_status crpx_alignment_find_taxon(const crpx_alignment *alignment, const char *name, size_t *taxon);
size_t crpx_alignment_patterns(const crpx_alignment *alignment);
/* Pattern data are column-major: data[pattern*taxa + taxon]. Missing lengths or
 * ragged data produce no pattern arrays; returned pointers remain valid until free. */
const unsigned char *crpx_alignment_pattern_data(const crpx_alignment *alignment);
const size_t *crpx_alignment_pattern_weights(const crpx_alignment *alignment);
const size_t *crpx_alignment_site_patterns(const crpx_alignment *alignment);
/* IUPAC DNA state mask A=1 C=2 G=4 T/U=8, ambiguity=union, -/?/X/N=15;
 * unrecognised characters return zero. */
uint8_t crpx_dna_state_mask(unsigned char symbol);
/* Pairwise deletion of all non-ACGT/U sites; counts original sites even when
 * patterns are compacted. Saturated JC/K2P returns +INFINITY with OK; no usable
 * sites returns NOT_FOUND. Invalid DNA symbols return INVALID_ARGUMENT. */
crpx_status crpx_alignment_pairwise_distance(const crpx_alignment *alignment, size_t a, size_t b,
                                            crpx_sequence_distance model, double *distance, size_t *used_sites);
/* The finite-distance matrix stores saturated estimates and pairs without
 * usable sites as NaN; use pairwise_distance to distinguish those cases. */
crpx_status crpx_alignment_distance_matrix(const crpx_alignment *alignment, crpx_sequence_distance model,
                                          crpx_distance_matrix_t *out);
/* Per-pattern 4-state tip vectors; caller allocates >= 4*patterns doubles.
 * Ambiguities have 1 for each allowed state. */
crpx_status crpx_alignment_tip_states(const crpx_alignment *alignment, size_t taxon, double *states, size_t capacity);

#ifdef __cplusplus
}
#endif
#endif
