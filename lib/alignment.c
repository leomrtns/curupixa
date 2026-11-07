/* Reimplementation of biomcmc-lib alignment (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "alignment.h"
#include "char_vector.h"
#include "file_compression.h"
#include "hashtable.h"
#include "nexus_common.h"
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct crpx_alignment {
  crpx_char_vector_t names, sequences;
  crpx_hashtable_t name_index;
  size_t taxa, sites, patterns;
  bool aligned;
  unsigned char *pattern_data;
  size_t *weights, *site_patterns;
};

static crpx_status make_patterns(crpx_alignment *a, bool compact)
{
  crpx_hashtable_t map = NULL;
  crpx_status status;
  unsigned char *column;
  if (!a->aligned || !a->sites) return CRPX_OK;
  if (a->taxa > SIZE_MAX / a->sites || a->sites > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  a->pattern_data = malloc(a->taxa * a->sites);
  a->weights = calloc(a->sites, sizeof(size_t));
  a->site_patterns = malloc(a->sites * sizeof(size_t));
  column = malloc(a->taxa);
  if (!a->pattern_data || !a->weights || !a->site_patterns || !column) { free(column); return CRPX_ERROR_ALLOCATION; }
  if (compact) {
    status = crpx_hashtable_new(a->sites, 0, &map);
    if (status != CRPX_OK) { free(column); return status; }
  }
  for (size_t site = 0; site < a->sites; ++site) {
    uint64_t pattern = a->patterns;
    for (size_t taxon = 0; taxon < a->taxa; ++taxon) column[taxon] = (unsigned char)crpx_alignment_sequence(a, taxon)[site];
    if (!compact || crpx_hashtable_get(map, column, a->taxa, &pattern) == CRPX_ERROR_NOT_FOUND) {
      pattern = a->patterns++;
      memcpy(a->pattern_data + (size_t)pattern * a->taxa, column, a->taxa);
      if (compact) {
        status = crpx_hashtable_put(map, column, a->taxa, pattern, NULL);
        if (status != CRPX_OK) { free(column); crpx_hashtable_free(map); return status; }
      }
    }
    a->site_patterns[site] = (size_t)pattern;
    ++a->weights[pattern];
  }
  free(column);
  crpx_hashtable_free(map);
  /* Release column capacity after pattern compression; failure to shrink is harmless. */
  unsigned char *data = realloc(a->pattern_data, a->taxa * a->patterns);
  if (data) a->pattern_data = data;
  size_t *weights = realloc(a->weights, a->patterns * sizeof(size_t));
  if (weights) a->weights = weights;
  return CRPX_OK;
}

crpx_status crpx_alignment_new(const char *const *names, const char *const *sequences, size_t taxa,
                               bool compact, crpx_alignment **out)
{
  crpx_alignment *a;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!names || !sequences || !taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  a = calloc(1, sizeof(*a));
  if (!a) return CRPX_ERROR_ALLOCATION;
  a->taxa = taxa;
  a->aligned = true;
  status = crpx_char_vector_new(taxa, &a->names);
  if (status == CRPX_OK) status = crpx_char_vector_new(taxa, &a->sequences);
  if (status == CRPX_OK) status = crpx_hashtable_new(taxa, 0, &a->name_index);
  for (size_t i = 0; i < taxa && status == CRPX_OK; ++i) {
    bool replaced;
    if (!names[i] || !*names[i] || !sequences[i]) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
    for (const unsigned char *p = (const unsigned char *)names[i]; *p; ++p)
      if (*p < 32 || *p == 127) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
    if (status != CRPX_OK) break;
    status = crpx_hashtable_put_string(a->name_index, names[i], i, &replaced);
    if (status != CRPX_OK) break;
    if (replaced) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
    status = crpx_char_vector_set(a->names, i, names[i]);
    if (status == CRPX_OK) status = crpx_char_vector_set(a->sequences, i, sequences[i]);
    if (status != CRPX_OK) break;
    /* This copy is privately owned until the immutable object is published. */
    char *sequence = (char *)crpx_char_vector_get(a->sequences, i);
    for (unsigned char *p = (unsigned char *)sequence; *p; ++p) {
      if (*p <= 32 || *p >= 127) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
      if (*p >= 'a' && *p <= 'z') *p = (unsigned char)(*p - 'a' + 'A');
    }
    size_t length = strlen(sequence);
    if (!i) a->sites = length;
    else if (length != a->sites) a->aligned = false;
  }
  if (status == CRPX_OK && compact && !a->aligned) status = CRPX_ERROR_INVALID_ARGUMENT;
  if (!a->aligned) a->sites = 0;
  if (status == CRPX_OK) status = make_patterns(a, compact);
  if (status != CRPX_OK) { crpx_alignment_free(a); return status; }
  *out = a;
  return CRPX_OK;
}

void crpx_alignment_free(crpx_alignment *a)
{
  if (!a) return;
  crpx_char_vector_free(a->names); crpx_char_vector_free(a->sequences); crpx_hashtable_free(a->name_index);
  free(a->pattern_data); free(a->weights); free(a->site_patterns); free(a);
}
size_t crpx_alignment_taxa(const crpx_alignment *a) { return a ? a->taxa : 0; }
bool crpx_alignment_is_aligned(const crpx_alignment *a) { return a && a->aligned; }
size_t crpx_alignment_sites(const crpx_alignment *a) { return a ? a->sites : 0; }
const char *crpx_alignment_name(const crpx_alignment *a, size_t taxon) { return a ? crpx_char_vector_get(a->names, taxon) : NULL; }
const char *crpx_alignment_sequence(const crpx_alignment *a, size_t taxon) { return a ? crpx_char_vector_get(a->sequences, taxon) : NULL; }
size_t crpx_alignment_patterns(const crpx_alignment *a) { return a ? a->patterns : 0; }
const unsigned char *crpx_alignment_pattern_data(const crpx_alignment *a) { return a ? a->pattern_data : NULL; }
const size_t *crpx_alignment_pattern_weights(const crpx_alignment *a) { return a ? a->weights : NULL; }
const size_t *crpx_alignment_site_patterns(const crpx_alignment *a) { return a ? a->site_patterns : NULL; }

crpx_status crpx_alignment_find_taxon(const crpx_alignment *a, const char *name, size_t *taxon)
{
  uint64_t value;
  if (!a || !taxon) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_hashtable_get_string(a->name_index, name, &value);
  if (status == CRPX_OK) *taxon = (size_t)value;
  return status;
}

uint8_t crpx_dna_state_mask(unsigned char symbol)
{
  if (symbol >= 'a' && symbol <= 'z') symbol = (unsigned char)(symbol - 'a' + 'A');
  switch (symbol) {
    case 'A': return 1; case 'C': return 2; case 'G': return 4; case 'T': case 'U': return 8;
    case 'R': return 5; case 'Y': return 10; case 'S': return 6; case 'W': return 9;
    case 'K': return 12; case 'M': return 3; case 'B': return 14; case 'D': return 13;
    case 'H': return 11; case 'V': return 7; case 'N': case 'X': case '?': case '-': return 15;
    default: return 0;
  }
}

crpx_status crpx_alignment_pairwise_distance(const crpx_alignment *a, size_t x, size_t y,
                                            crpx_sequence_distance model, double *distance, size_t *used)
{
  size_t valid = 0, transitions = 0, transversions = 0;
  if (!a || !distance || !used || !a->aligned || model < CRPX_DISTANCE_P || model > CRPX_DISTANCE_K2P)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (x >= a->taxa || y >= a->taxa) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t p = 0; p < a->patterns; ++p) {
    uint8_t left = crpx_dna_state_mask(a->pattern_data[p*a->taxa+x]);
    uint8_t right = crpx_dna_state_mask(a->pattern_data[p*a->taxa+y]);
    if (!left || !right) return CRPX_ERROR_INVALID_ARGUMENT;
    if ((left & (left - 1)) || (right & (right - 1))) continue;
    size_t weight = a->weights[p];
    valid += weight;
    if (left != right) {
      if ((left | right) == 5 || (left | right) == 10) transitions += weight;
      else transversions += weight;
    }
  }
  if (!valid) return CRPX_ERROR_NOT_FOUND;
  double p = (double)transitions / (double)valid, q = (double)transversions / (double)valid;
  double result = p + q;
  if (model == CRPX_DISTANCE_JC69) result = result >= 0.75 ? INFINITY : -0.75 * log1p(-4.0 * result / 3.0);
  if (model == CRPX_DISTANCE_K2P)
    result = 2*p+q >= 1 || 2*q >= 1 ? INFINITY : -0.5*log1p(-2*p-q) - 0.25*log1p(-2*q);
  *distance = result;
  *used = valid;
  return CRPX_OK;
}

crpx_status crpx_alignment_distance_matrix(const crpx_alignment *a, crpx_sequence_distance model,
                                          crpx_distance_matrix_t *out)
{
  crpx_distance_matrix_t matrix;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!a || !a->aligned || model < CRPX_DISTANCE_P || model > CRPX_DISTANCE_K2P) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_distance_matrix_new(a->taxa, &matrix);
  if (status != CRPX_OK) return status;
  for (size_t i = 1; i < a->taxa && status == CRPX_OK; ++i) for (size_t j = 0; j < i; ++j) {
    double distance;
    size_t used;
    status = crpx_alignment_pairwise_distance(a, i, j, model, &distance, &used);
    if (status == CRPX_ERROR_NOT_FOUND) { status = CRPX_OK; distance = NAN; }
    if (status != CRPX_OK) break;
    /* Saturated model distances are retained as missing in a finite-distance
     * matrix; use pairwise_distance to distinguish saturation from no sites. */
    status = crpx_distance_matrix_set(matrix, i, j, isfinite(distance) ? distance : NAN);
    if (status != CRPX_OK) break;
  }
  if (status != CRPX_OK) { crpx_distance_matrix_free(matrix); return status; }
  *out = matrix;
  return CRPX_OK;
}

crpx_status crpx_alignment_tip_states(const crpx_alignment *a, size_t taxon, double *states, size_t capacity)
{
  if (!a || !a->aligned || (!states && a->patterns)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxon >= a->taxa) return CRPX_ERROR_OUT_OF_RANGE;
  if (a->patterns > SIZE_MAX / 4) return CRPX_ERROR_OVERFLOW;
  if (capacity < a->patterns * 4) return CRPX_ERROR_OUT_OF_RANGE;
  for (size_t p = 0; p < a->patterns; ++p)
    if (!crpx_dna_state_mask(a->pattern_data[p*a->taxa+taxon])) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t p = 0; p < a->patterns; ++p) {
    uint8_t mask = crpx_dna_state_mask(a->pattern_data[p*a->taxa+taxon]);
    for (size_t state = 0; state < 4; ++state) states[p*4+state] = (mask & (1u << state)) ? 1 : 0;
  }
  return CRPX_OK;
}

static crpx_status vectors_to_alignment(crpx_char_vector_t names, crpx_char_vector_t seq, bool compact,
                                       crpx_alignment **out)
{
  size_t count = crpx_char_vector_size(names);
  if (!count) return CRPX_ERROR_PARSE;
  if (count > SIZE_MAX / sizeof(char *)) return CRPX_ERROR_OVERFLOW;
  const char **n = malloc(count * sizeof(*n)), **s = malloc(count * sizeof(*s));
  if (!n || !s) { free(n); free(s); return CRPX_ERROR_ALLOCATION; }
  for (size_t i = 0; i < count; ++i) { n[i] = crpx_char_vector_get(names, i); s[i] = crpx_char_vector_get(seq, i); }
  crpx_status status = crpx_alignment_new(n, s, count, compact, out);
  free(n); free(s);
  return status;
}

static crpx_status parse_fasta(char *text, bool compact, crpx_alignment **out)
{
  crpx_char_vector_t names = NULL, sequences = NULL;
  crpx_status status = crpx_char_vector_new(0, &names);
  if (status == CRPX_OK) status = crpx_char_vector_new(0, &sequences);
  for (char *line = text; *line && status == CRPX_OK;) {
    char *next = strchr(line, '\n');
    if (next) *next++ = '\0'; else next = line + strlen(line);
    while (isspace((unsigned char)*line)) ++line;
    size_t length = strlen(line);
    while (length && isspace((unsigned char)line[length-1])) line[--length] = '\0';
    if (*line == '>') {
      ++line; while (isspace((unsigned char)*line)) ++line;
      if (!*line) { status = CRPX_ERROR_PARSE; break; }
      status = crpx_char_vector_push(names, line);
      if (status == CRPX_OK) status = crpx_char_vector_push(sequences, "");
    } else if (*line && *line != ';' && *line != '#') {
      if (!crpx_char_vector_size(names)) { status = CRPX_ERROR_PARSE; break; }
      char *to = line;
      for (char *from = line; *from; ++from) if (!isspace((unsigned char)*from)) *to++ = *from;
      *to = '\0';
      status = crpx_char_vector_append(sequences, crpx_char_vector_size(sequences)-1, line);
    }
    line = next;
  }
  if (status == CRPX_OK) status = vectors_to_alignment(names, sequences, compact, out);
  crpx_char_vector_free(names); crpx_char_vector_free(sequences);
  return status;
}

/* Bounded tokeniser for NEXUS headers and matrix names; edits only its private
 * input. Newline handling for sequential/interleaved rows stays in parse_nexus. */
static char *next_token(char **cursor, crpx_status *status)
{
  char *p = *cursor;
  while (*p && isspace((unsigned char)*p)) ++p;
  if (!*p) { *cursor = p; return NULL; }
  char *start = p;
  if (*p == '\'') {
    char *to = start;
    bool closed = false;
    ++p;
    while (*p) {
      if (*p == '\'') {
        if (p[1] == '\'') { *to++ = '\''; p += 2; }
        else { ++p; closed = true; break; }
      } else *to++ = *p++;
    }
    if (!closed || (*p && !isspace((unsigned char)*p))) { *status = CRPX_ERROR_PARSE; return NULL; }
    *to = '\0';
    *cursor = p;
    return start;
  }
  while (*p && !isspace((unsigned char)*p)) ++p;
  if (*p) *p++ = '\0';
  *cursor = p;
  return start;
}

static int equal_word(const char *a, const char *b)
{
  for (; *a && *b; ++a, ++b) if (toupper((unsigned char)*a) != toupper((unsigned char)*b)) return 0;
  return !*a && !*b;
}

/* Tokenise a command before MATRIX without modifying delimiters. */
static crpx_status header_token(const char **cursor, char word[64])
{
  const char *p = *cursor;
  while (isspace((unsigned char)*p)) p++;
  size_t length = 0;
  if (*p == '=') { word[length++] = *p++; }
  else if (*p == '\'') {
    p++;
    bool closed = false;
    while (*p) {
      if (*p == '\'') {
        if (p[1] == '\'') p++;
        else { p++; closed = true; break; }
      }
      if (length == 63) return CRPX_ERROR_UNSUPPORTED;
      word[length++] = *p++;
    }
    if (!closed || (*p && !isspace((unsigned char)*p))) return CRPX_ERROR_PARSE;
  } else {
    while (*p && !isspace((unsigned char)*p) && *p != '=') {
      if (length == 63) return CRPX_ERROR_UNSUPPORTED;
      word[length++] = *p++;
    }
  }
  word[length] = 0;
  *cursor = p;
  return length ? CRPX_OK : CRPX_ERROR_NOT_FOUND;
}

static crpx_status parse_dimensions(const char *cursor, size_t *taxa, size_t *sites)
{
  char key[64], equals[64], number[64];
  crpx_status status;
  while ((status = header_token(&cursor, key)) == CRPX_OK) {
    if (equal_word(key, "NEWTAXA")) continue;
    if (!equal_word(key, "NTAX") && !equal_word(key, "NCHAR")) return CRPX_ERROR_UNSUPPORTED;
    if (header_token(&cursor, equals) != CRPX_OK || strcmp(equals, "=") ||
        header_token(&cursor, number) != CRPX_OK) return CRPX_ERROR_PARSE;
    size_t value = 0;
    for (size_t i = 0; number[i]; i++) {
      if (number[i] < '0' || number[i] > '9') return CRPX_ERROR_PARSE;
      unsigned digit = (unsigned)(number[i] - '0');
      if (value > (SIZE_MAX - digit) / 10) return CRPX_ERROR_OVERFLOW;
      value = value * 10 + digit;
    }
    if (!value) return CRPX_ERROR_PARSE;
    size_t *target = equal_word(key, "NTAX") ? taxa : sites;
    if (*target && *target != value) return CRPX_ERROR_PARSE;
    *target = value;
  }
  return status == CRPX_ERROR_NOT_FOUND ? CRPX_OK : status;
}

static crpx_status parse_format(const char *cursor, bool *interleaved)
{
  char key[64], equals[64], value[64];
  crpx_status status;
  while ((status = header_token(&cursor, key)) == CRPX_OK) {
    bool layout_flag = equal_word(key, "INTERLEAVE") || equal_word(key, "LABELS");
    if (!layout_flag && !equal_word(key, "DATATYPE") && !equal_word(key, "GAP") &&
        !equal_word(key, "MISSING") && !equal_word(key, "MATCHCHAR")) return CRPX_ERROR_UNSUPPORTED;
    const char *before = cursor;
    status = header_token(&cursor, equals);
    if (status != CRPX_OK && status != CRPX_ERROR_NOT_FOUND) return status;
    if (layout_flag && (status == CRPX_ERROR_NOT_FOUND || strcmp(equals, "="))) {
      cursor = before;
      if (equal_word(key, "INTERLEAVE")) *interleaved = true;
      continue;
    }
    if (status != CRPX_OK || strcmp(equals, "=") || header_token(&cursor, value) != CRPX_OK)
      return CRPX_ERROR_PARSE;
    if (layout_flag) {
      if (!equal_word(value, "YES") && !equal_word(value, "NO")) return CRPX_ERROR_PARSE;
      if (equal_word(key, "LABELS") && equal_word(value, "NO")) return CRPX_ERROR_UNSUPPORTED;
      if (equal_word(key, "INTERLEAVE")) *interleaved = equal_word(value, "YES");
    } else if (equal_word(key, "DATATYPE")) {
      if (!equal_word(value, "DNA") && !equal_word(value, "RNA") && !equal_word(value, "NUCLEOTIDE"))
        return CRPX_ERROR_UNSUPPORTED;
    } else {
      char expected = equal_word(key, "GAP") ? '-' : equal_word(key, "MISSING") ? '?' : '.';
      if (value[0] != expected || value[1]) return CRPX_ERROR_UNSUPPORTED;
    }
  }
  return status == CRPX_ERROR_NOT_FOUND ? CRPX_OK : status;
}

static crpx_status parse_nexus(const char *text, bool compact, crpx_alignment **out)
{
  char *copy = NULL, *matrix = NULL;
  crpx_char_vector_t names = NULL, sequences = NULL;
  crpx_status status = crpx_nexus_strip_comments(text, &copy, NULL);
  size_t expected_taxa = 0, expected_sites = 0;
  bool interleaved = false, alignment_closed = false;
  int block = 0; /* DATA/CHARACTERS=1, TAXA=2, other=3, outside=0 */
  if (status != CRPX_OK) return status;
  char *cursor = copy + 6; /* from_text already checked the #NEXUS header */
  while (*cursor) {
    while (isspace((unsigned char)*cursor)) cursor++;
    if (!*cursor) break;
    char *end = cursor;
    bool quoted = false;
    for (; *end; end++) {
      if (*end == '\'') { if (quoted && end[1] == '\'') end++; else quoted = !quoted; }
      else if (!quoted && *end == ';') break;
    }
    if (!*end || quoted) { status = CRPX_ERROR_PARSE; goto done; }
    *end = 0;
    const char *arguments = cursor;
    char command[64];
    status = header_token(&arguments, command);
    if (status != CRPX_OK) { if (status == CRPX_ERROR_NOT_FOUND) status = CRPX_ERROR_PARSE; goto done; }
    if (equal_word(command, "BEGIN")) {
      char name[64], extra[64];
      if (block || header_token(&arguments, name) != CRPX_OK || header_token(&arguments, extra) != CRPX_ERROR_NOT_FOUND) {
        status = CRPX_ERROR_PARSE; goto done;
      }
      block = equal_word(name, "DATA") || equal_word(name, "CHARACTERS") ? 1 : equal_word(name, "TAXA") ? 2 : 3;
    } else if (equal_word(command, "END") || equal_word(command, "ENDBLOCK")) {
      if (!block) { status = CRPX_ERROR_PARSE; goto done; }
      if (block == 1 && matrix) alignment_closed = true;
      block = 0;
    } else if (block == 1 || block == 2) {
      if (block == 1 && matrix) { status = CRPX_ERROR_UNSUPPORTED; goto done; }
      if (equal_word(command, "DIMENSIONS")) status = parse_dimensions(arguments, &expected_taxa, &expected_sites);
      else if (block == 1 && equal_word(command, "FORMAT")) status = parse_format(arguments, &interleaved);
      else if (block == 1 && equal_word(command, "MATRIX")) {
        if (matrix) { status = CRPX_ERROR_UNSUPPORTED; goto done; }
        if (!expected_taxa || !expected_sites) { status = CRPX_ERROR_PARSE; goto done; }
        matrix = (char *)arguments;
      } else if (!equal_word(command, "TITLE") && !equal_word(command, "LINK") &&
                 !(block == 2 && equal_word(command, "TAXLABELS"))) status = CRPX_ERROR_UNSUPPORTED;
      if (status != CRPX_OK) goto done;
    } else if (!block) { status = CRPX_ERROR_PARSE; goto done; }
    cursor = end + 1;
  }
  if (!matrix || !alignment_closed || block) { status = CRPX_ERROR_PARSE; goto done; }
  status = crpx_char_vector_new(0, &names);
  if (status == CRPX_OK) status = crpx_char_vector_new(0, &sequences);
  size_t current = SIZE_MAX;
  for (char *line = matrix; *line && status == CRPX_OK;) {
    char *next = strchr(line, '\n');
    if (next) *next++ = '\0'; else next = line + strlen(line);
    char *part = line;
    while (isspace((unsigned char)*part)) part++;
    bool quoted_label = *part == '\'';
    char *name = next_token(&part, &status);
    if (status != CRPX_OK) break;
    if (!name) { line = next; continue; }
    char *sequence = part;
    while (isspace((unsigned char)*sequence)) ++sequence;
    bool continuation = !interleaved && current != SIZE_MAX &&
                        crpx_char_vector_length(sequences, current) < expected_sites;
    size_t named = SIZE_MAX;
    if (continuation && crpx_char_vector_find(names, name, &named) == CRPX_OK) continuation = false;
    if (!*sequence || continuation) {
      /* A continuation of the preceding sequential row. */
      if (interleaved) { status = CRPX_ERROR_UNSUPPORTED; break; }
      if (quoted_label) { status = CRPX_ERROR_PARSE; break; }
      if (current == SIZE_MAX || crpx_char_vector_length(sequences, current) >= expected_sites) {
        status = CRPX_ERROR_UNSUPPORTED; break;
      }
      if (continuation && *sequence) {
        size_t prefix = strlen(name);
        name[prefix] = ' '; /* rejoin the token separator for whitespace removal below */
      }
      sequence = name;
    } else {
      status = crpx_char_vector_find(names, name, &current);
      if (status == CRPX_ERROR_NOT_FOUND) {
        current = crpx_char_vector_size(names);
        if (current >= expected_taxa) { status = CRPX_ERROR_PARSE; break; }
        status = crpx_char_vector_push(names, name);
        if (status == CRPX_OK) status = crpx_char_vector_push(sequences, "");
      }
      if (status != CRPX_OK) break;
    }
    char *to = sequence;
    for (char *p = sequence; *p; ++p) if (!isspace((unsigned char)*p)) *to++ = *p;
    *to = '\0';
    status = crpx_char_vector_append(sequences, current, sequence);
    if (status == CRPX_OK && crpx_char_vector_length(sequences, current) > expected_sites) status = CRPX_ERROR_PARSE;
    line = next;
  }
  if (status == CRPX_OK && crpx_char_vector_size(names) != expected_taxa) status = CRPX_ERROR_PARSE;
  for (size_t i = 0; status == CRPX_OK && i < expected_taxa; ++i) {
    if (crpx_char_vector_length(sequences, i) != expected_sites) { status = CRPX_ERROR_PARSE; break; }
    char *sequence = (char *)crpx_char_vector_get(sequences, i);
    const char *first = crpx_char_vector_get(sequences, 0);
    for (size_t j = 0; j < expected_sites; ++j) if (sequence[j] == '.') {
      if (!i) { status = CRPX_ERROR_PARSE; break; }
      sequence[j] = first[j];
    }
  }
  if (status == CRPX_OK) status = vectors_to_alignment(names, sequences, compact, out);
done:
  free(copy); crpx_char_vector_free(names); crpx_char_vector_free(sequences);
  return status;
}

crpx_status crpx_alignment_from_text(const char *text, bool compact, crpx_alignment **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!text) return CRPX_ERROR_INVALID_ARGUMENT;
  while (isspace((unsigned char)*text)) ++text;
  if (*text == '#') {
    char tag[7] = {0};
    size_t i = 0;
    while (i < 6 && text[i]) { tag[i] = (char)toupper((unsigned char)text[i]); ++i; }
    if (!strcmp(tag, "#NEXUS") && (!text[6] || isspace((unsigned char)text[6])))
      return parse_nexus(text, compact, out);
  }
  /* FASTA permits leading comment lines as well as comments between records. */
  while (*text == ';' || *text == '#') {
    const char *next = strchr(text, '\n');
    if (!next) return CRPX_ERROR_PARSE;
    text = next + 1;
    while (isspace((unsigned char)*text)) ++text;
  }
  if (*text == '>') {
    size_t length = strlen(text);
    if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
    char *copy = malloc(length + 1);
    if (!copy) return CRPX_ERROR_ALLOCATION;
    memcpy(copy, text, length + 1);
    crpx_status status = parse_fasta(copy, compact, out);
    free(copy); return status;
  }
  return CRPX_ERROR_PARSE;
}

crpx_status crpx_alignment_read(const char *path, bool compact, crpx_alignment **out)
{
  crpx_file *file;
  crpx_char_vector_t text;
  char *line = NULL;
  size_t capacity = 0, length;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  crpx_status status = crpx_file_open(path, "r", CRPX_FILE_AUTO, &file);
  if (status != CRPX_OK) return status;
  status = crpx_char_vector_new(1, &text);
  while (status == CRPX_OK) {
    status = crpx_file_getline(file, &line, &capacity, &length);
    if (status != CRPX_OK) break;
    if (memchr(line, 0, length)) { status = CRPX_ERROR_PARSE; break; }
    status = crpx_char_vector_append(text, 0, line);
  }
  if (status == CRPX_ERROR_NOT_FOUND) status = CRPX_OK;
  crpx_status closed = crpx_file_close(file);
  if (status == CRPX_OK) status = closed;
  if (status == CRPX_OK) status = crpx_alignment_from_text(crpx_char_vector_get(text, 0), compact, out);
  free(line); crpx_char_vector_free(text);
  return status;
}

crpx_status crpx_alignment_write_fasta(const crpx_alignment *a, const char *path, size_t width)
{
  crpx_file *file;
  if (!a || !width) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_file_open(path, "w", CRPX_FILE_AUTO, &file);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < a->taxa && status == CRPX_OK; ++i) {
    const char *name = crpx_alignment_name(a, i), *sequence = crpx_alignment_sequence(a, i);
    status = crpx_file_write(file, ">", 1);
    if (status == CRPX_OK) status = crpx_file_write(file, name, strlen(name));
    if (status == CRPX_OK) status = crpx_file_write(file, "\n", 1);
    size_t length = strlen(sequence);
    for (size_t position = 0; position < length && status == CRPX_OK;) {
      size_t count = length-position < width ? length-position : width;
      status = crpx_file_write(file, sequence+position, count);
      if (status == CRPX_OK) status = crpx_file_write(file, "\n", 1);
      position += count;
    }
  }
  crpx_status closed = crpx_file_close(file);
  return status == CRPX_OK ? closed : status;
}
