/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Format reference: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md */
#include "gff3_format.h"
#include "file_compression.h"
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct attribute { char *key, **values; size_t count; };
struct feature { crpx_gff3_feature view; struct attribute *attributes; size_t count; };
struct record { char *raw; size_t feature; };
struct crpx_gff3 {
  struct feature *features;
  size_t feature_count;
  crpx_gff3_region *regions;
  size_t region_count;
  struct record *records;
  size_t record_count;
  crpx_alignment *sequences;
  char **fasta_ids;
};
struct buffer { char *data; size_t length, capacity; };

static crpx_status append(struct buffer *buffer, const char *text, size_t length)
{
  if (length > SIZE_MAX - buffer->length - 1) return CRPX_ERROR_OVERFLOW;
  size_t required = buffer->length + length + 1;
  if (required > buffer->capacity) {
    size_t capacity = buffer->capacity ? buffer->capacity : 256;
    while (capacity < required) capacity = capacity > SIZE_MAX / 2 ? required : capacity * 2;
    char *replacement = realloc(buffer->data, capacity);
    if (!replacement) return CRPX_ERROR_ALLOCATION;
    buffer->data = replacement; buffer->capacity = capacity;
  }
  if (length) memcpy(buffer->data + buffer->length, text, length);
  buffer->length += length; buffer->data[buffer->length] = 0;
  return CRPX_OK;
}

static crpx_status grow(void **array, size_t count, size_t item)
{
  if (count && (count & (count - 1))) return CRPX_OK;
  if (count > SIZE_MAX / 2) return CRPX_ERROR_OVERFLOW;
  size_t capacity = count ? count * 2 : 1;
  if (capacity > SIZE_MAX / item) return CRPX_ERROR_OVERFLOW;
  void *replacement = realloc(*array, capacity * item);
  if (!replacement) return CRPX_ERROR_ALLOCATION;
  *array = replacement;
  return CRPX_OK;
}

static int hex_value(unsigned char character)
{
  if (character >= '0' && character <= '9') return character - '0';
  if (character >= 'a' && character <= 'f') return character - 'a' + 10;
  if (character >= 'A' && character <= 'F') return character - 'A' + 10;
  return -1;
}

static bool seqid_character(unsigned char character)
{
  return (character >= 'a' && character <= 'z') || (character >= 'A' && character <= 'Z') ||
         (character >= '0' && character <= '9') || strchr(".:^*$@!+_?-|", character);
}

static bool valid_raw_seqid(const char *text)
{
  for (size_t i = 0; text[i]; i++) {
    if (text[i] == '%') {
      if (!text[i + 1] || !text[i + 2]) return false;
      i += 2;
    } else if (!seqid_character((unsigned char)text[i])) return false;
  }
  return true;
}

static crpx_status decode(const char *text, size_t length, char **result)
{
  if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  char *value = malloc(length + 1);
  if (!value) return CRPX_ERROR_ALLOCATION;
  size_t written = 0;
  for (size_t i = 0; i < length; i++) {
    unsigned char character = (unsigned char)text[i];
    if (character == '%') {
      if (length - i < 3 || hex_value((unsigned char)text[i + 1]) < 0 || hex_value((unsigned char)text[i + 2]) < 0) {
        free(value); return CRPX_ERROR_PARSE;
      }
      character = (unsigned char)((hex_value((unsigned char)text[i + 1]) << 4) | hex_value((unsigned char)text[i + 2]));
      i += 2;
      if (!character) { free(value); return CRPX_ERROR_UNSUPPORTED; }
    } else if (character < 32 || character == 127) { free(value); return CRPX_ERROR_PARSE; }
    value[written++] = (char)character;
  }
  value[written] = 0; *result = value;
  return CRPX_OK;
}

static void free_feature(struct feature *feature)
{
  free((char *)feature->view.seqid); free((char *)feature->view.source); free((char *)feature->view.type);
  for (size_t i = 0; i < feature->count; i++) {
    free(feature->attributes[i].key);
    for (size_t j = 0; j < feature->attributes[i].count; j++) free(feature->attributes[i].values[j]);
    free(feature->attributes[i].values);
  }
  free(feature->attributes);
}

void crpx_gff3_free(crpx_gff3 *document)
{
  if (!document) return;
  for (size_t i = 0; i < document->feature_count; i++) free_feature(&document->features[i]);
  for (size_t i = 0; i < document->region_count; i++) free((char *)document->regions[i].seqid);
  for (size_t i = 0; i < document->record_count; i++) free(document->records[i].raw);
  if (document->fasta_ids) for (size_t i = 0; i < crpx_alignment_taxa(document->sequences); i++) free(document->fasta_ids[i]);
  free(document->fasta_ids); crpx_alignment_free(document->sequences);
  free(document->features); free(document->regions); free(document->records); free(document);
}

static crpx_status integer(const char *text, uint64_t *result)
{
  if (!*text) return CRPX_ERROR_PARSE;
  uint64_t value = 0;
  for (; *text; text++) {
    if (*text < '0' || *text > '9') return CRPX_ERROR_PARSE;
    unsigned digit = (unsigned)(*text - '0');
    if (value > (UINT64_MAX - digit) / 10) return CRPX_ERROR_OVERFLOW;
    value = value * 10 + digit;
  }
  if (!value) return CRPX_ERROR_PARSE;
  *result = value; return CRPX_OK;
}

static crpx_status attributes(struct feature *feature, char *text)
{
  if (!strcmp(text, ".")) return CRPX_OK;
  if (!*text) return CRPX_ERROR_PARSE;
  for (char *cursor = text; *cursor;) {
    char *end = strchr(cursor, ';');
    if (!end) end = cursor + strlen(cursor);
    char *equals = memchr(cursor, '=', (size_t)(end - cursor));
    if (!equals || equals == cursor) return CRPX_ERROR_PARSE;
    crpx_status status = grow((void **)&feature->attributes, feature->count, sizeof(*feature->attributes));
    if (status != CRPX_OK) return status;
    struct attribute *attribute = &feature->attributes[feature->count++];
    memset(attribute, 0, sizeof(*attribute));
    status = decode(cursor, (size_t)(equals - cursor), &attribute->key);
    if (status != CRPX_OK) return status;
    for (size_t i = 0; i + 1 < feature->count; i++)
      if (!strcmp(feature->attributes[i].key, attribute->key)) return CRPX_ERROR_PARSE;
    char *value = equals + 1;
    for (;;) {
      char *comma = memchr(value, ',', (size_t)(end - value));
      char *stop = comma ? comma : end;
      status = grow((void **)&attribute->values, attribute->count, sizeof(*attribute->values));
      if (status != CRPX_OK) return status;
      for (char *raw = value; raw < stop; raw++) if (*raw == '=' || *raw == '&') return CRPX_ERROR_PARSE;
      char *decoded = NULL;
      status = decode(value, (size_t)(stop - value), &decoded);
      if (status != CRPX_OK) return status;
      attribute->values[attribute->count++] = decoded;
      if (!comma) break;
      value = comma + 1;
    }
    if ((!strcmp(attribute->key, "ID") || !strcmp(attribute->key, "Parent")))
      for (size_t i = 0; i < attribute->count; i++) if (!*attribute->values[i]) return CRPX_ERROR_PARSE;
    cursor = *end ? end + 1 : end;
  }
  return CRPX_OK;
}

static bool decimal_digit(char character)
{
  return character >= '0' && character <= '9';
}

/* Decimal/scientific syntax only. Conversion and formatting never change or inspect
 * shared locale state. Keep 19 significant digits, as the Newick reader does. */
static crpx_status parse_score(const char *text, double *result)
{
  bool negative = false, decimal = false, any_digit = false, significant = false, exponent_negative = false;
  unsigned kept = 0;
  long double mantissa = 0;
  int64_t fractional = 0, skipped = 0, exponent = 0;
  if (*text == '+' || *text == '-') negative = *text++ == '-';
  for (;;) {
    char character = *text;
    if (decimal_digit(character)) {
      any_digit = true;
      if (decimal) fractional++;
      if (character != '0') significant = true;
      if (significant) {
        if (kept < 19) { mantissa = mantissa * 10 + (character - '0'); kept++; }
        else skipped++;
      }
      text++;
    } else if (character == '.' && !decimal) { decimal = true; text++; }
    else break;
  }
  if (!any_digit) return CRPX_ERROR_PARSE;
  if (*text == 'e' || *text == 'E') {
    text++;
    if (*text == '+' || *text == '-') exponent_negative = *text++ == '-';
    if (!decimal_digit(*text)) return CRPX_ERROR_PARSE;
    while (decimal_digit(*text)) {
      if (exponent < 1000000) exponent = exponent * 10 + (*text - '0');
      text++;
    }
  }
  if (*text) return CRPX_ERROR_PARSE;
  if (exponent_negative) exponent = -exponent;
  exponent += skipped - fractional;
  long double value = mantissa;
  /* Scale in chunks: a lone 10^-324 factor could underflow prematurely on systems
   * where long double has the same range as double, although 5e-324 is representable. */
  while (exponent > 100 && isfinite(value) && value != 0) { value *= 1e100L; exponent -= 100; }
  while (exponent < -100 && value != 0) { value *= 1e-100L; exponent += 100; }
  if (value != 0 && isfinite(value)) value *= powl(10.0L, (long double)exponent);
  if (negative) value = -value;
  double converted = (double)value;
  if (!isfinite(converted) || (mantissa != 0 && converted == 0)) return CRPX_ERROR_PARSE;
  *result = converted;
  return CRPX_OK;
}

static crpx_status format_score(double value, char result[64])
{
  char number[64];
  int count = snprintf(number, sizeof(number), "%.17g", value);
  if (!isfinite(value) || count < 0 || (size_t)count >= sizeof(number)) return CRPX_ERROR_OVERFLOW;
  size_t source = 0, target = 0;
  while (number[source]) {
    if (decimal_digit(number[source]) || number[source] == '+' || number[source] == '-' ||
        number[source] == 'e' || number[source] == 'E') result[target++] = number[source++];
    else {
      result[target++] = '.';
      do { source++; }
      while (number[source] && !decimal_digit(number[source]) && number[source] != 'e' && number[source] != 'E');
    }
  }
  result[target] = 0;
  return CRPX_OK;
}

static crpx_status parse_feature(char *line, struct feature *feature)
{
  char *fields[9]; fields[0] = line;
  size_t count = 1;
  for (char *cursor = line; *cursor; cursor++) if (*cursor == '\t') {
    if (count == 9) return CRPX_ERROR_PARSE;
    *cursor = 0; fields[count++] = cursor + 1;
  }
  if (count != 9) return CRPX_ERROR_PARSE;
  for (size_t i = 0; i < 9; i++) if (!*fields[i]) return CRPX_ERROR_PARSE;
  if (!valid_raw_seqid(fields[0])) return CRPX_ERROR_PARSE;
  crpx_status status = decode(fields[0], strlen(fields[0]), (char **)&feature->view.seqid);
  if (status == CRPX_OK) status = decode(fields[1], strlen(fields[1]), (char **)&feature->view.source);
  if (status == CRPX_OK) status = decode(fields[2], strlen(fields[2]), (char **)&feature->view.type);
  if (status == CRPX_OK) status = integer(fields[3], &feature->view.start);
  if (status == CRPX_OK) status = integer(fields[4], &feature->view.end);
  if (status != CRPX_OK) return status;
  if (feature->view.start > feature->view.end || !strcmp(fields[0], ".") || !strcmp(fields[2], "."))
    return CRPX_ERROR_PARSE;
  feature->view.start--;
  if (strcmp(fields[5], ".")) {
    status = parse_score(fields[5], &feature->view.score);
    if (status != CRPX_OK) return status;
    feature->view.has_score = true;
  }
  if (fields[6][1] || !strchr("+-.?", fields[6][0])) return CRPX_ERROR_PARSE;
  feature->view.strand = fields[6][0];
  feature->view.phase = -1;
  if (strcmp(fields[7], ".")) {
    if (fields[7][1] || fields[7][0] < '0' || fields[7][0] > '2') return CRPX_ERROR_PARSE;
    feature->view.phase = fields[7][0] - '0';
  }
  if (!strcmp(feature->view.type, "CDS") && feature->view.phase == -1) return CRPX_ERROR_PARSE;
  return attributes(feature, fields[8]);
}

static char *token(char **cursor)
{
  char *start = *cursor;
  while (*start == ' ' || *start == '\t') start++;
  char *end = start;
  while (*end && *end != ' ' && *end != '\t') end++;
  if (*end) *end++ = 0;
  *cursor = end;
  return *start ? start : NULL;
}

static crpx_status parse_region(crpx_gff3 *document, const char *line)
{
  char *copy = malloc(strlen(line) + 1);
  if (!copy) return CRPX_ERROR_ALLOCATION;
  strcpy(copy, line);
  char *cursor = copy + 17;
  char *name = token(&cursor), *start = token(&cursor), *end = token(&cursor);
  crpx_gff3_region region = {0};
  crpx_status status = !name || !start || !end || token(&cursor) ? CRPX_ERROR_PARSE : CRPX_OK;
  if (status == CRPX_OK && !valid_raw_seqid(name)) status = CRPX_ERROR_PARSE;
  if (status == CRPX_OK) status = decode(name, strlen(name), (char **)&region.seqid);
  if (status == CRPX_OK) status = integer(start, &region.start);
  if (status == CRPX_OK) status = integer(end, &region.end);
  if (status == CRPX_OK && region.start > region.end) status = CRPX_ERROR_PARSE;
  if (status == CRPX_OK) {
    region.start--;
    for (size_t i = 0; i < document->region_count; i++)
      if (!strcmp(document->regions[i].seqid, region.seqid)) status = CRPX_ERROR_PARSE;
  }
  if (status == CRPX_OK) status = grow((void **)&document->regions, document->region_count, sizeof(region));
  if (status == CRPX_OK) document->regions[document->region_count++] = region;
  else free((char *)region.seqid);
  free(copy); return status;
}

crpx_status crpx_gff3_set_sequences(crpx_gff3 *document, const crpx_alignment *sequences)
{
  if (!document || !sequences) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t count = crpx_alignment_taxa(sequences);
  if (count > SIZE_MAX / sizeof(char *)) return CRPX_ERROR_OVERFLOW;
  const char **names = calloc(count ? count : 1, sizeof(*names));
  const char **strings = calloc(count ? count : 1, sizeof(*strings));
  char **identifiers = calloc(count ? count : 1, sizeof(*identifiers));
  if (!names || !strings || !identifiers) { free(names); free(strings); free(identifiers); return CRPX_ERROR_ALLOCATION; }
  crpx_status status = CRPX_OK;
  for (size_t i = 0; i < count && status == CRPX_OK; i++) {
    names[i] = crpx_alignment_name(sequences, i);
    strings[i] = crpx_alignment_sequence(sequences, i);
    size_t length = 0;
    while (names[i][length] && !isspace((unsigned char)names[i][length])) length++;
    status = decode(names[i], length, &identifiers[i]);
    if (status == CRPX_OK && !length) status = CRPX_ERROR_PARSE;
    for (size_t j = 0; j < i && status == CRPX_OK; j++)
      if (!strcmp(identifiers[i], identifiers[j])) status = CRPX_ERROR_PARSE;
  }
  crpx_alignment *copy = NULL;
  if (status == CRPX_OK) status = crpx_alignment_new(names, strings, count, false, &copy);
  free(names); free(strings);
  if (status != CRPX_OK) {
    for (size_t i = 0; i < count; i++) free(identifiers[i]);
    free(identifiers); return status;
  }
  if (document->fasta_ids) for (size_t i = 0; i < crpx_alignment_taxa(document->sequences); i++) free(document->fasta_ids[i]);
  free(document->fasta_ids); crpx_alignment_free(document->sequences);
  document->sequences = copy; document->fasta_ids = identifiers;
  return CRPX_OK;
}

crpx_status crpx_gff3_from_text(const char *text, crpx_gff3 **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!text) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_gff3 *document = calloc(1, sizeof(*document));
  if (!document) return CRPX_ERROR_ALLOCATION;
  crpx_status status = CRPX_OK;
  bool version = false;
  while (*text && status == CRPX_OK) {
    const char *end = strchr(text, '\n');
    if (!end) end = text + strlen(text);
    size_t length = (size_t)(end - text);
    if (length && text[length - 1] == '\r') length--;
    char *line = malloc(length + 1);
    if (!line) { status = CRPX_ERROR_ALLOCATION; break; }
    memcpy(line, text, length); line[length] = 0;
    text = *end ? end + 1 : end;
    if (!version && *line) {
      if (strcmp(line, "##gff-version 3")) { free(line); status = CRPX_ERROR_PARSE; break; }
      version = true;
    }
    if (!strncmp(line, "##gff-version", 13) && strcmp(line, "##gff-version 3")) {
      free(line); status = CRPX_ERROR_UNSUPPORTED; break;
    }
    if (!strcmp(line, "##FASTA")) {
      free(line);
      crpx_alignment *alignment = NULL;
      status = crpx_alignment_from_text(text, false, &alignment);
      if (status == CRPX_OK) status = crpx_gff3_set_sequences(document, alignment);
      crpx_alignment_free(alignment);
      break;
    }
    status = grow((void **)&document->records, document->record_count, sizeof(*document->records));
    if (status != CRPX_OK) { free(line); break; }
    struct record record = {line, 0};
    if (*line && *line != '#') {
      struct feature feature = {0};
      status = parse_feature(line, &feature);
      if (status == CRPX_OK) status = grow((void **)&document->features, document->feature_count, sizeof(feature));
      if (status == CRPX_OK) {
        record.feature = document->feature_count;
        document->features[document->feature_count++] = feature;
      } else free_feature(&feature);
      free(line); record.raw = NULL;
    } else if (!strncmp(line, "##sequence-region", 17) && (!line[17] || isspace((unsigned char)line[17])))
      status = parse_region(document, line);
    if (status == CRPX_OK) document->records[document->record_count++] = record;
    else free(record.raw);
  }
  if (status == CRPX_OK && !version) status = CRPX_ERROR_PARSE;
  if (status != CRPX_OK) { crpx_gff3_free(document); return status; }
  *result = document; return CRPX_OK;
}

size_t crpx_gff3_feature_count(const crpx_gff3 *document) { return document ? document->feature_count : 0; }
size_t crpx_gff3_region_count(const crpx_gff3 *document) { return document ? document->region_count : 0; }
const crpx_gff3_feature *crpx_gff3_feature_at(const crpx_gff3 *document, size_t feature)
{ return document && feature < document->feature_count ? &document->features[feature].view : NULL; }
const crpx_gff3_region *crpx_gff3_region_at(const crpx_gff3 *document, size_t region)
{ return document && region < document->region_count ? &document->regions[region] : NULL; }
size_t crpx_gff3_attribute_count(const crpx_gff3 *document, size_t feature)
{ return document && feature < document->feature_count ? document->features[feature].count : 0; }
static const struct attribute *attribute_at(const crpx_gff3 *document, size_t feature, size_t attribute)
{
  if (!document || feature >= document->feature_count || attribute >= document->features[feature].count) return NULL;
  return &document->features[feature].attributes[attribute];
}
const char *crpx_gff3_attribute_key(const crpx_gff3 *document, size_t feature, size_t attribute)
{ const struct attribute *item = attribute_at(document, feature, attribute); return item ? item->key : NULL; }
size_t crpx_gff3_attribute_value_count(const crpx_gff3 *document, size_t feature, size_t attribute)
{ const struct attribute *item = attribute_at(document, feature, attribute); return item ? item->count : 0; }
const char *crpx_gff3_attribute_value(const crpx_gff3 *document, size_t feature, size_t attribute, size_t value)
{ const struct attribute *item = attribute_at(document, feature, attribute); return item && value < item->count ? item->values[value] : NULL; }
crpx_status crpx_gff3_find_attribute(const crpx_gff3 *document, size_t feature, const char *key, size_t *attribute)
{
  if (!document || feature >= document->feature_count || !key || !attribute) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < document->features[feature].count; i++) if (!strcmp(document->features[feature].attributes[i].key, key)) {
    *attribute = i; return CRPX_OK;
  }
  return CRPX_ERROR_NOT_FOUND;
}
const crpx_alignment *crpx_gff3_sequences(const crpx_gff3 *document) { return document ? document->sequences : NULL; }
crpx_status crpx_gff3_sequence(const crpx_gff3 *document, const char *seqid, const char **sequence, size_t *length)
{
  if (!document || !seqid || !sequence || !length) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < crpx_alignment_taxa(document->sequences); i++) if (!strcmp(document->fasta_ids[i], seqid)) {
    *sequence = crpx_alignment_sequence(document->sequences, i); *length = strlen(*sequence); return CRPX_OK;
  }
  return CRPX_ERROR_NOT_FOUND;
}
crpx_status crpx_gff3_query(const crpx_gff3 *document, const char *seqid, uint64_t start, uint64_t end,
                          const char *type, size_t *indices, size_t capacity, size_t *count)
{
  if (!document || !count || start >= end || (!indices && capacity)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t matches = 0;
  for (size_t i = 0; i < document->feature_count; i++) {
    const crpx_gff3_feature *feature = &document->features[i].view;
    if (feature->start < end && start < feature->end && (!seqid || !strcmp(feature->seqid, seqid)) &&
        (!type || !strcmp(feature->type, type))) matches++;
  }
  *count = matches;
  if (capacity < matches) return CRPX_ERROR_OUT_OF_RANGE;
  matches = 0;
  for (size_t i = 0; i < document->feature_count; i++) {
    const crpx_gff3_feature *feature = &document->features[i].view;
    if (feature->start < end && start < feature->end && (!seqid || !strcmp(feature->seqid, seqid)) &&
        (!type || !strcmp(feature->type, type))) indices[matches++] = i;
  }
  return CRPX_OK;
}

static crpx_status escape(struct buffer *buffer, const char *text, bool seqid, bool attribute)
{
  static const char hex[] = "0123456789ABCDEF";
  for (const unsigned char *cursor = (const unsigned char *)text; *cursor; cursor++) {
    unsigned char value = *cursor;
    bool reserved = value < 32 || value == 127 || value == '%' ||
                    (attribute && strchr(";=,&", value)) ||
                    (seqid && !seqid_character(value));
    char encoded[3] = {'%', hex[value >> 4], hex[value & 15]};
    crpx_status status = append(buffer, reserved ? encoded : (const char *)cursor, reserved ? 3 : 1);
    if (status != CRPX_OK) return status;
  }
  return CRPX_OK;
}

crpx_status crpx_gff3_to_text(const crpx_gff3 *document, char **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!document) return CRPX_ERROR_INVALID_ARGUMENT;
  struct buffer output = {0};
  crpx_status status = CRPX_OK;
  for (size_t i = 0; i < document->record_count && status == CRPX_OK; i++) {
    if (document->records[i].raw) status = append(&output, document->records[i].raw, strlen(document->records[i].raw));
    else {
      const struct feature *feature = &document->features[document->records[i].feature];
      const crpx_gff3_feature *view = &feature->view;
      status = escape(&output, view->seqid, true, false);
      if (status == CRPX_OK) status = append(&output, "\t", 1);
      if (status == CRPX_OK) status = escape(&output, view->source, false, false);
      if (status == CRPX_OK) status = append(&output, "\t", 1);
      if (status == CRPX_OK) status = escape(&output, view->type, false, false);
      char values[160], score[64], phase[2] = {'.', 0};
      if (view->has_score) {
        if (status == CRPX_OK) status = format_score(view->score, score);
        if (status != CRPX_OK) break;
      } else strcpy(score, ".");
      if (view->phase >= 0) phase[0] = (char)('0' + view->phase);
      int length = snprintf(values, sizeof(values), "\t%" PRIu64 "\t%" PRIu64 "\t%s\t%c\t%s\t",
                            view->start + 1, view->end, score, view->strand, phase);
      if (status == CRPX_OK) status = append(&output, values, (size_t)length);
      if (!feature->count && status == CRPX_OK) status = append(&output, ".", 1);
      for (size_t j = 0; j < feature->count && status == CRPX_OK; j++) {
        if (j) status = append(&output, ";", 1);
        if (status == CRPX_OK) status = escape(&output, feature->attributes[j].key, false, true);
        if (status == CRPX_OK) status = append(&output, "=", 1);
        for (size_t k = 0; k < feature->attributes[j].count && status == CRPX_OK; k++) {
          if (k) status = append(&output, ",", 1);
          if (status == CRPX_OK) status = escape(&output, feature->attributes[j].values[k], false, true);
        }
      }
    }
    if (status == CRPX_OK) status = append(&output, "\n", 1);
  }
  if (document->sequences && status == CRPX_OK) status = append(&output, "##FASTA\n", 8);
  for (size_t i = 0; i < crpx_alignment_taxa(document->sequences) && status == CRPX_OK; i++) {
    const char *name = crpx_alignment_name(document->sequences, i), *sequence = crpx_alignment_sequence(document->sequences, i);
    status = append(&output, ">", 1);
    if (status == CRPX_OK) status = append(&output, name, strlen(name));
    if (status == CRPX_OK) status = append(&output, "\n", 1);
    if (status == CRPX_OK) status = append(&output, sequence, strlen(sequence));
    if (status == CRPX_OK) status = append(&output, "\n", 1);
  }
  if (status != CRPX_OK) { free(output.data); return status; }
  *result = output.data; return CRPX_OK;
}

crpx_status crpx_gff3_read(const char *path, crpx_gff3 **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!path) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_file *file;
  crpx_status status = crpx_file_open(path, "r", CRPX_FILE_AUTO, &file);
  if (status != CRPX_OK) return status;
  struct buffer input = {0};
  char block[8192];
  while (status == CRPX_OK) {
    size_t count;
    status = crpx_file_read(file, block, sizeof(block), &count);
    if (status != CRPX_OK || !count) break;
    if (memchr(block, 0, count)) { status = CRPX_ERROR_PARSE; break; }
    status = append(&input, block, count);
  }
  crpx_status closed = crpx_file_close(file);
  if (status == CRPX_OK) status = closed;
  if (status == CRPX_OK) status = crpx_gff3_from_text(input.data ? input.data : "", result);
  free(input.data); return status;
}

crpx_status crpx_gff3_write(const crpx_gff3 *document, const char *path)
{
  if (!document || !path) return CRPX_ERROR_INVALID_ARGUMENT;
  char *text;
  crpx_status status = crpx_gff3_to_text(document, &text);
  if (status != CRPX_OK) return status;
  crpx_file *file;
  status = crpx_file_open(path, "w", CRPX_FILE_AUTO, &file);
  if (status == CRPX_OK) {
    status = crpx_file_write(file, text, strlen(text));
    crpx_status closed = crpx_file_close(file);
    if (status == CRPX_OK) status = closed;
  }
  free(text); return status;
}
