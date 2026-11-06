/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2022-today Leonardo de Oliveira Martins
 * Reimplementation of biomcmc-lib Newick/tree collections with explicit ownership.
 */
#include "newick_space.h"
#include "file_compression.h"
#include "topology_distance.h"
#include "hashtable.h"
#include <stdlib.h>
#include <string.h>

typedef struct { crpx_tree *tree; uint64_t weight; } space_entry;

struct crpx_tree_space {
  space_entry *entries;
  size_t size, capacity;
  uint64_t total_weight;
};

crpx_status
crpx_tree_space_new(crpx_tree_space **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = (crpx_tree_space *) calloc(1, sizeof(**result));
  return *result ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}

void
crpx_tree_space_free(crpx_tree_space *space)
{
  size_t i;
  if (!space) return;
  for (i = 0; i < space->size; ++i) crpx_tree_free(space->entries[i].tree);
  free(space->entries);
  free(space);
}

size_t crpx_tree_space_size(const crpx_tree_space *space) { return space ? space->size : 0; }
uint64_t crpx_tree_space_total_weight(const crpx_tree_space *space) { return space ? space->total_weight : 0; }
const crpx_tree *crpx_tree_space_get(const crpx_tree_space *space, size_t index)
{ return space && index < space->size ? space->entries[index].tree : NULL; }

crpx_status
crpx_tree_space_weight(const crpx_tree_space *space, size_t index, uint64_t *weight)
{
  if (!space || !weight) return CRPX_ERROR_INVALID_ARGUMENT;
  if (index >= space->size) return CRPX_ERROR_OUT_OF_RANGE;
  *weight = space->entries[index].weight;
  return CRPX_OK;
}

static crpx_status
space_reserve(crpx_tree_space *space, size_t size)
{
  size_t capacity;
  space_entry *replacement;
  if (size <= space->capacity) return CRPX_OK;
  capacity = space->capacity ? space->capacity : 16;
  while (capacity < size) {
    if (capacity > SIZE_MAX / 2) { capacity = size; break; }
    capacity *= 2;
  }
  if (capacity > SIZE_MAX / sizeof(*replacement)) return CRPX_ERROR_OVERFLOW;
  replacement = (space_entry *) realloc(space->entries, capacity * sizeof(*replacement));
  if (!replacement) return CRPX_ERROR_ALLOCATION;
  space->entries = replacement;
  space->capacity = capacity;
  return CRPX_OK;
}

static crpx_status
space_take(crpx_tree_space *space, crpx_tree *tree, uint64_t weight)
{
  crpx_status status;
  if (!weight) return CRPX_ERROR_INVALID_ARGUMENT;
  if (space->size == SIZE_MAX || UINT64_MAX - space->total_weight < weight) return CRPX_ERROR_OVERFLOW;
  status = space_reserve(space, space->size + 1);
  if (status != CRPX_OK) return status;
  space->entries[space->size].tree = tree;
  space->entries[space->size++].weight = weight;
  space->total_weight += weight;
  return CRPX_OK;
}

crpx_status
crpx_tree_space_add(crpx_tree_space *space, const crpx_tree *tree, uint64_t weight)
{
  crpx_tree *copy;
  crpx_status status;
  if (!space || !tree || !weight) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!crpx_tree_node_count(tree)) return CRPX_ERROR_TOPOLOGY;
  if (UINT64_MAX - space->total_weight < weight) return CRPX_ERROR_OVERFLOW;
  status = crpx_tree_clone(tree, &copy);
  if (status != CRPX_OK) return status;
  status = space_take(space, copy, weight);
  if (status != CRPX_OK) crpx_tree_free(copy);
  return status;
}

crpx_status
crpx_tree_space_clone(const crpx_tree_space *space, crpx_tree_space **result)
{
  crpx_tree_space *copy;
  crpx_status status;
  size_t i;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!space) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_space_new(&copy);
  if (status != CRPX_OK) return status;
  for (i = 0; i < space->size; ++i) {
    status = crpx_tree_space_add(copy, space->entries[i].tree, space->entries[i].weight);
    if (status != CRPX_OK) { crpx_tree_space_free(copy); return status; }
  }
  *result = copy;
  return CRPX_OK;
}

typedef struct {
  crpx_tree_space *space;
  char *text;
  size_t size, capacity, offset, start, comments, error_offset;
  bool quoted, meaningful;
  unsigned format;
  bool nexus_header, in_block, in_trees, saw_trees;
  crpx_hashtable_t translation;
  char **translated_labels;
  size_t translated_count, translated_capacity;
} newick_stream;

static bool
stream_space(char c)
{
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

static bool
nexus_equal(const char *first, const char *second)
{
  unsigned char a, b;
  if (!first || !second) return false;
  do {
    a = (unsigned char) *first++; b = (unsigned char) *second++;
    if (a >= 'A' && a <= 'Z') a = (unsigned char) (a + ('a' - 'A'));
    if (b >= 'A' && b <= 'Z') b = (unsigned char) (b + ('a' - 'A'));
    if (a != b) return false;
  } while (a);
  return true;
}

static crpx_status
nexus_skip(const char *text, size_t *position)
{
  size_t depth;
  for (;;) {
    while (stream_space(text[*position])) ++*position;
    if (text[*position] != '[') return CRPX_OK;
    depth = 1; ++*position;
    while (depth) {
      if (!text[*position]) return CRPX_ERROR_PARSE;
      if (text[*position] == '[') ++depth;
      else if (text[*position] == ']') --depth;
      ++*position;
    }
  }
}

static bool
nexus_punctuation(char c)
{
  return c == ',' || c == ';' || c == '=' || c == '*';
}

static crpx_status
nexus_token(const char *text, size_t *position, char **result)
{
  size_t start, end, length = 0, i;
  char *token;
  crpx_status status;
  *result = NULL;
  status = nexus_skip(text, position);
  if (status != CRPX_OK || !text[*position]) return status;
  if (text[*position] == '\'') {
    start = ++*position;
    end = start;
    for (;;) {
      if (!text[end]) { *position = end; return CRPX_ERROR_PARSE; }
      if (text[end] == '\'') {
        ++end;
        if (text[end] != '\'') break;
      }
      ++length; ++end;
    }
    token = (char *) malloc(length + 1);
    if (!token) return CRPX_ERROR_ALLOCATION;
    i = 0;
    while (start + 1 < end) {
      token[i++] = text[start++];
      if (token[i - 1] == '\'') ++start;
    }
    token[i] = '\0';
    *position = end;
  }
  else {
    start = *position;
    if (nexus_punctuation(text[*position])) ++*position;
    else while (text[*position] && !stream_space(text[*position]) && text[*position] != '[' &&
                !nexus_punctuation(text[*position])) ++*position;
    length = *position - start;
    token = (char *) malloc(length + 1);
    if (!token) return CRPX_ERROR_ALLOCATION;
    memcpy(token, text + start, length); token[length] = '\0';
  }
  *result = token;
  return CRPX_OK;
}

static void
clear_translation(newick_stream *stream)
{
  size_t i;
  crpx_hashtable_free(stream->translation); stream->translation = NULL;
  for (i = 0; i < stream->translated_count; ++i) free(stream->translated_labels[i]);
  free(stream->translated_labels); stream->translated_labels = NULL;
  stream->translated_count = stream->translated_capacity = 0;
}

static crpx_status
parse_translation(newick_stream *stream, size_t *position)
{
  char *key = NULL, *value = NULL, *separator = NULL, **replacement;
  size_t capacity;
  bool replaced;
  crpx_status status;
  if (stream->translation) return CRPX_ERROR_PARSE;
  status = crpx_hashtable_new(16, 0, &stream->translation);
  if (status != CRPX_OK) return status;
  for (;;) {
    status = nexus_token(stream->text, position, &key);
    if (status != CRPX_OK) break;
    if (!key || !*key || nexus_punctuation(key[0])) { status = CRPX_ERROR_PARSE; break; }
    status = nexus_skip(stream->text, position);
    if (status != CRPX_OK) break;
    if (stream->text[*position] != '\'' && nexus_punctuation(stream->text[*position])) {
      status = CRPX_ERROR_PARSE; break;
    }
    status = nexus_token(stream->text, position, &value);
    if (status != CRPX_OK) break;
    if (!value) { status = CRPX_ERROR_PARSE; break; }
    if (stream->translated_count == stream->translated_capacity) {
      capacity = stream->translated_capacity ? stream->translated_capacity * 2 : 16;
      if (capacity < stream->translated_capacity || capacity > SIZE_MAX / sizeof(*replacement)) {
        status = CRPX_ERROR_OVERFLOW; break;
      }
      replacement = (char **) realloc(stream->translated_labels, capacity * sizeof(*replacement));
      if (!replacement) { status = CRPX_ERROR_ALLOCATION; break; }
      stream->translated_labels = replacement; stream->translated_capacity = capacity;
    }
    status = crpx_hashtable_put_string(stream->translation, key, stream->translated_count, &replaced);
    if (status != CRPX_OK) break;
    if (replaced) { status = CRPX_ERROR_PARSE; break; }
    stream->translated_labels[stream->translated_count++] = value; value = NULL;
    free(key); key = NULL;
    status = nexus_token(stream->text, position, &separator);
    if (status != CRPX_OK) break;
    if (nexus_equal(separator, ";")) break;
    if (!nexus_equal(separator, ",")) { status = CRPX_ERROR_PARSE; break; }
    free(separator); separator = NULL;
  }
  free(key); free(value); free(separator);
  return status;
}

static crpx_status
parse_nexus_statement(newick_stream *stream)
{
  char *command = NULL, *token = NULL;
  size_t position = 0, local_offset, i;
  uint64_t translation;
  const char *label;
  crpx_tree *tree = NULL;
  crpx_status status = nexus_token(stream->text, &position, &command);
  if (status != CRPX_OK) goto finished;
  if (!stream->nexus_header) {
    if (!nexus_equal(command, "#nexus")) { status = CRPX_ERROR_PARSE; goto finished; }
    stream->nexus_header = true;
    free(command); command = NULL;
    status = nexus_token(stream->text, &position, &command);
    if (status != CRPX_OK) goto finished;
  }
  if (nexus_equal(command, "begin")) {
    if (stream->in_block) { status = CRPX_ERROR_PARSE; goto finished; }
    status = nexus_token(stream->text, &position, &token);
    if (status != CRPX_OK) goto finished;
    if (!token || nexus_equal(token, ";")) { status = CRPX_ERROR_PARSE; goto finished; }
    stream->in_block = true;
    stream->in_trees = nexus_equal(token, "trees");
    if (stream->in_trees) { clear_translation(stream); stream->saw_trees = true; }
    free(token); token = NULL;
    status = nexus_token(stream->text, &position, &token);
    if (status == CRPX_OK && !nexus_equal(token, ";")) status = CRPX_ERROR_PARSE;
  }
  else if (nexus_equal(command, "end") || nexus_equal(command, "endblock")) {
    if (!stream->in_block) { status = CRPX_ERROR_PARSE; goto finished; }
    stream->in_block = stream->in_trees = false;
    status = nexus_token(stream->text, &position, &token);
    if (status == CRPX_OK && !nexus_equal(token, ";")) status = CRPX_ERROR_PARSE;
  }
  else if (stream->in_trees && nexus_equal(command, "translate")) status = parse_translation(stream, &position);
  else if (stream->in_trees && (nexus_equal(command, "tree") || nexus_equal(command, "utree"))) {
    status = nexus_token(stream->text, &position, &token);
    if (status != CRPX_OK) goto finished;
    if (nexus_equal(token, "*")) {
      free(token); token = NULL;
      status = nexus_token(stream->text, &position, &token);
      if (status != CRPX_OK) goto finished;
    }
    if (!token || nexus_equal(token, "=") || nexus_equal(token, ";")) { status = CRPX_ERROR_PARSE; goto finished; }
    free(token); token = NULL;
    status = nexus_token(stream->text, &position, &token);
    if (status != CRPX_OK) goto finished;
    if (!nexus_equal(token, "=")) { status = CRPX_ERROR_PARSE; goto finished; }
    status = crpx_tree_from_newick(stream->text + position, &tree, &local_offset);
    if (status != CRPX_OK) { position += local_offset; goto finished; }
    if (stream->translation) for (i = 0; i < crpx_tree_node_count(tree); ++i) {
      if (crpx_tree_first_child(tree, (crpx_node_id) i) != CRPX_NODE_NONE) continue;
      label = crpx_tree_label(tree, (crpx_node_id) i);
      if (!label) { status = CRPX_ERROR_PARSE; goto finished; }
      status = crpx_hashtable_get_string(stream->translation, label, &translation);
      if (status == CRPX_ERROR_NOT_FOUND) status = CRPX_ERROR_PARSE;
      if (status != CRPX_OK) goto finished;
      status = crpx_tree_set_label(tree, (crpx_node_id) i, stream->translated_labels[translation]);
      if (status != CRPX_OK) goto finished;
    }
    status = space_take(stream->space, tree, 1);
    if (status == CRPX_OK) tree = NULL;
  }
finished:
  free(command); free(token); crpx_tree_free(tree);
  if (status != CRPX_OK) stream->error_offset = stream->start + position;
  return status;
}

static crpx_status
stream_feed(newick_stream *stream, const char *text, size_t length)
{
  size_t i, capacity, local_offset;
  char c, *replacement, *first_token;
  crpx_tree *tree;
  crpx_status status;
  for (i = 0; i < length; ++i) {
    c = text[i];
    if (!c) { stream->error_offset = stream->offset; return CRPX_ERROR_PARSE; }
    if (stream->size > SIZE_MAX - 2 || stream->offset == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
    if (stream->size + 2 > stream->capacity) {
      capacity = stream->capacity ? stream->capacity : 256;
      if (capacity <= stream->size + 1) {
        if (capacity > SIZE_MAX / 2) return CRPX_ERROR_OVERFLOW;
        capacity *= 2;
      }
      replacement = (char *) realloc(stream->text, capacity);
      if (!replacement) return CRPX_ERROR_ALLOCATION;
      stream->text = replacement; stream->capacity = capacity;
    }
    stream->text[stream->size++] = c;
    stream->text[stream->size] = '\0';
    ++stream->offset;
    if (stream->comments) {
      if (c == '[') ++stream->comments;
      else if (c == ']') --stream->comments;
      continue;
    }
    if (c == '\'') { stream->quoted = !stream->quoted; stream->meaningful = true; continue; }
    if (stream->quoted) continue;
    if (c == '[') { stream->comments = 1; continue; }
    /* A legal Newick label may itself begin with '#'. Inspect the complete first token
     * before selecting NEXUS, while retaining streaming operation across header chunks.
     */
    if (!stream->format && !stream_space(c)) stream->format = c == '#' ? 3 : 1;
    if (c != ';') { if (!stream_space(c)) stream->meaningful = true; continue; }
    if (stream->format == 3) {
      local_offset = 0;
      status = nexus_token(stream->text, &local_offset, &first_token);
      if (status != CRPX_OK) return status;
      stream->format = nexus_equal(first_token, "#nexus") ? 2 : 1;
      free(first_token);
    }
    if (stream->format == 2) {
      status = parse_nexus_statement(stream);
      if (status != CRPX_OK) return status;
    }
    else {
      status = crpx_tree_from_newick(stream->text, &tree, &local_offset);
      if (status != CRPX_OK) { stream->error_offset = stream->start + local_offset; return status; }
      status = space_take(stream->space, tree, 1);
      if (status != CRPX_OK) { crpx_tree_free(tree); return status; }
    }
    stream->size = 0; stream->start = stream->offset; stream->meaningful = false;
  }
  stream->error_offset = stream->offset;
  return CRPX_OK;
}

static crpx_status
stream_finish(const newick_stream *stream)
{
  return stream->quoted || stream->comments || stream->meaningful ||
         (stream->format == 2 && (!stream->nexus_header || !stream->saw_trees || stream->in_block)) ?
         CRPX_ERROR_PARSE : CRPX_OK;
}

crpx_status
crpx_tree_space_add_trees(crpx_tree_space *space, const char *text, size_t *error_offset)
{
  crpx_tree_space *incoming;
  newick_stream stream = { 0 };
  crpx_status status;
  if (error_offset) *error_offset = 0;
  if (!space || !text) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_space_new(&incoming);
  if (status != CRPX_OK) return status;
  stream.space = incoming;
  status = stream_feed(&stream, text, strlen(text));
  if (status == CRPX_OK) status = stream_finish(&stream);
  if (error_offset) *error_offset = stream.error_offset;
  free(stream.text);
  clear_translation(&stream);
  if (status == CRPX_OK) {
    if (SIZE_MAX - space->size < incoming->size || UINT64_MAX - space->total_weight < incoming->total_weight)
      status = CRPX_ERROR_OVERFLOW;
    else status = space_reserve(space, space->size + incoming->size);
    if (status == CRPX_OK && incoming->size) {
      memcpy(space->entries + space->size, incoming->entries, incoming->size * sizeof(*space->entries));
      space->size += incoming->size;
      space->total_weight += incoming->total_weight;
      incoming->size = 0;
    }
  }
  crpx_tree_space_free(incoming);
  return status;
}

crpx_status
crpx_tree_space_add_newick(crpx_tree_space *space, const char *text, size_t *error_offset)
{ return crpx_tree_space_add_trees(space, text, error_offset); }

crpx_status
crpx_tree_space_read_file(const char *path, crpx_tree_space **result, size_t *error_offset)
{
  crpx_tree_space *space;
  crpx_file *file = NULL;
  crpx_status status, close_status;
  newick_stream stream = { 0 };
  char buffer[4096];
  size_t count;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (error_offset) *error_offset = 0;
  if (!path) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_space_new(&space);
  if (status != CRPX_OK) return status;
  stream.space = space;
  status = crpx_file_open(path, "r", CRPX_FILE_AUTO, &file);
  if (status == CRPX_OK) for (;;) {
    status = crpx_file_read(file, buffer, sizeof(buffer), &count);
    if (status != CRPX_OK) break;
    if (!count) { status = stream_finish(&stream); break; }
    status = stream_feed(&stream, buffer, count);
    if (status != CRPX_OK) break;
  }
  if (file) {
    close_status = crpx_file_close(file);
    if (status == CRPX_OK) status = close_status;
  }
  if (error_offset) *error_offset = stream.error_offset;
  free(stream.text);
  clear_translation(&stream);
  if (status != CRPX_OK) crpx_tree_space_free(space);
  else *result = space;
  return status;
}

crpx_status
crpx_tree_space_write_file(const crpx_tree_space *space, const char *path)
{
  crpx_file *file;
  char *text;
  size_t i;
  crpx_status status, close_status;
  if (!space || !path) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_file_open(path, "w", CRPX_FILE_AUTO, &file);
  if (status != CRPX_OK) return status;
  for (i = 0; i < space->size; ++i) {
    status = crpx_tree_to_newick(space->entries[i].tree, &text);
    if (status != CRPX_OK) break;
    status = crpx_file_write(file, text, strlen(text));
    crpx_tree_string_free(text);
    if (status == CRPX_OK) status = crpx_file_write(file, "\n", 1);
    if (status != CRPX_OK) break;
  }
  close_status = crpx_file_close(file);
  return status == CRPX_OK ? close_status : status;
}

crpx_status
crpx_tree_space_compact(crpx_tree_space *space, bool rooted)
{
  crpx_tree_space *compacted;
  crpx_status status;
  size_t i, j, distance;
  space_entry *old_entries;
  size_t old_size, old_capacity;
  uint64_t old_weight;
  if (!space) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_space_new(&compacted);
  if (status != CRPX_OK) return status;
  for (i = 0; i < space->size; ++i) {
    for (j = 0; j < compacted->size; ++j) {
      status = crpx_tree_rf_distance(space->entries[i].tree, compacted->entries[j].tree, rooted, &distance);
      if (status != CRPX_OK) goto finished;
      if (!distance) break;
    }
    if (j == compacted->size) {
      /* Self-comparison also validates labels when adding the first representative. */
      if (!j) {
        status = crpx_tree_rf_distance(space->entries[i].tree, space->entries[i].tree, rooted, &distance);
        if (status != CRPX_OK) goto finished;
      }
      status = crpx_tree_space_add(compacted, space->entries[i].tree, space->entries[i].weight);
      if (status != CRPX_OK) goto finished;
    }
    else {
      compacted->entries[j].weight += space->entries[i].weight;
      compacted->total_weight += space->entries[i].weight;
    }
  }
  old_entries = space->entries; old_size = space->size; old_capacity = space->capacity; old_weight = space->total_weight;
  *space = *compacted;
  compacted->entries = old_entries; compacted->size = old_size;
  compacted->capacity = old_capacity; compacted->total_weight = old_weight;
finished:
  crpx_tree_space_free(compacted);
  return status;
}

crpx_status
crpx_tree_space_select(const crpx_tree_space *space, size_t burnin, size_t thin, crpx_tree_space **result)
{
  crpx_tree_space *selected;
  size_t i;
  crpx_status status;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!space || !thin) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_tree_space_new(&selected);
  if (status != CRPX_OK) return status;
  for (i = burnin; i < space->size;) {
    status = crpx_tree_space_add(selected, space->entries[i].tree, space->entries[i].weight);
    if (status != CRPX_OK) { crpx_tree_space_free(selected); return status; }
    if (SIZE_MAX - i < thin) break;
    i += thin;
  }
  *result = selected;
  return CRPX_OK;
}

crpx_status
crpx_tree_space_merge(crpx_tree_space *space, const crpx_tree_space *from, uint64_t multiplier)
{
  crpx_tree_space *incoming;
  crpx_status status;
  size_t i;
  if (!space || !from || !multiplier) return CRPX_ERROR_INVALID_ARGUMENT;
  if (from->total_weight > (UINT64_MAX - space->total_weight) / multiplier ||
      SIZE_MAX - space->size < from->size) return CRPX_ERROR_OVERFLOW;
  status = crpx_tree_space_clone(from, &incoming);
  if (status != CRPX_OK) return status;
  status = space_reserve(space, space->size + incoming->size);
  if (status == CRPX_OK && incoming->size) {
    for (i = 0; i < incoming->size; ++i) incoming->entries[i].weight *= multiplier;
    memcpy(space->entries + space->size, incoming->entries, incoming->size * sizeof(*space->entries));
    space->size += incoming->size;
    space->total_weight += incoming->total_weight * multiplier;
    incoming->size = 0;
  }
  crpx_tree_space_free(incoming);
  return status;
}

crpx_status
crpx_tree_space_frequency(const crpx_tree_space *space, size_t index, double *frequency)
{
  if (!space || !frequency) return CRPX_ERROR_INVALID_ARGUMENT;
  if (index >= space->size) return CRPX_ERROR_OUT_OF_RANGE;
  *frequency = (double) space->entries[index].weight / (double) space->total_weight;
  return CRPX_OK;
}
