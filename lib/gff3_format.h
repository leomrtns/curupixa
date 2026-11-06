/* Reimplementation of biomcmc-lib GFF3 support (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_GFF3_FORMAT_H
#define CURUPIXA_GFF3_FORMAT_H
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "alignment.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_gff3 crpx_gff3;
typedef struct {
  const char *seqid, *source, *type;
  uint64_t start, end; /* zero-based half-open interval; GFF start-1 and end */
  double score;
  bool has_score;
  char strand; /* +, -, . or ? */
  int phase; /* -1 for missing, otherwise 0..2 */
} crpx_gff3_feature;
typedef struct {
  const char *seqid;
  uint64_t start, end; /* zero-based half-open */
} crpx_gff3_region;

/* Owns parsed data. Requires ##gff-version 3. Records preserve input order, comments
 * and unknown directives. Decodes percent escapes; comma-separated attribute values
 * remain separate, so an escaped comma remains part of one value. %00 is unsupported.
 * Scores accept finite decimal/scientific notation and always use a decimal point,
 * independently of the calling process locale. This parser checks syntax; it does
 * not validate the Sequence Ontology or Parent graph. */
crpx_status crpx_gff3_from_text(const char *text, crpx_gff3 **result);
crpx_status crpx_gff3_read(const char *path, crpx_gff3 **result);
void crpx_gff3_free(crpx_gff3 *document);
crpx_status crpx_gff3_to_text(const crpx_gff3 *document, char **result);
crpx_status crpx_gff3_write(const crpx_gff3 *document, const char *path);
size_t crpx_gff3_feature_count(const crpx_gff3 *document);
const crpx_gff3_feature *crpx_gff3_feature_at(const crpx_gff3 *document, size_t index);
size_t crpx_gff3_region_count(const crpx_gff3 *document);
const crpx_gff3_region *crpx_gff3_region_at(const crpx_gff3 *document, size_t index);
size_t crpx_gff3_attribute_count(const crpx_gff3 *document, size_t feature);
const char *crpx_gff3_attribute_key(const crpx_gff3 *document, size_t feature, size_t attribute);
size_t crpx_gff3_attribute_value_count(const crpx_gff3 *document, size_t feature, size_t attribute);
const char *crpx_gff3_attribute_value(const crpx_gff3 *document, size_t feature, size_t attribute, size_t value);
crpx_status crpx_gff3_find_attribute(const crpx_gff3 *document, size_t feature, const char *key, size_t *attribute);
/* Borrowed optional FASTA, with descriptions retained in alignment names. */
const crpx_alignment *crpx_gff3_sequences(const crpx_gff3 *document);
/* Deep-copies sequences, linking exact first whitespace-delimited FASTA identifiers to
 * decoded seqids. Duplicate identifiers fail. Mutation requires exclusive ownership. */
crpx_status crpx_gff3_set_sequences(crpx_gff3 *document, const crpx_alignment *sequences);
crpx_status crpx_gff3_sequence(const crpx_gff3 *document, const char *seqid, const char **sequence, size_t *length);
/* Find all overlapping features, including nested/overlapping annotations. NULL seqid/type
 * matches all. count reports required capacity; OUT_OF_RANGE leaves indices untouched.
 * Immutable queries may run concurrently. Returned feature/string pointers last until free. */
crpx_status crpx_gff3_query(const crpx_gff3 *document, const char *seqid, uint64_t start, uint64_t end,
                          const char *type, size_t *indices, size_t capacity, size_t *count);
#ifdef __cplusplus
}
#endif
#endif
