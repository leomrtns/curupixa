/* Reimplementation of biomcmc-lib char_vector (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "char_vector.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct string_slot { char *text; size_t length, capacity; };
struct crpx_char_vector_struct { size_t size, capacity; struct string_slot *slots; };

static crpx_status reserve_slots(crpx_char_vector_t v, size_t size)
{
  struct string_slot *slots;
  if (size <= v->capacity) return CRPX_OK;
  if (size > SIZE_MAX / sizeof(*slots)) return CRPX_ERROR_OVERFLOW;
  size_t capacity = v->capacity <= SIZE_MAX / 2 ? v->capacity * 2 : size;
  if (capacity < size || capacity > SIZE_MAX / sizeof(*slots)) capacity = size;
  slots = realloc(v->slots, capacity * sizeof(*slots));
  if (!slots) return CRPX_ERROR_ALLOCATION;
  memset(slots + v->capacity, 0, (capacity - v->capacity) * sizeof(*slots));
  v->slots = slots;
  v->capacity = capacity;
  return CRPX_OK;
}

crpx_status crpx_char_vector_new(size_t size, crpx_char_vector_t *out)
{
  crpx_char_vector_t v;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  v = calloc(1, sizeof(*v));
  if (!v) return CRPX_ERROR_ALLOCATION;
  status = reserve_slots(v, size);
  if (status != CRPX_OK) { free(v); return status; }
  v->size = size;
  *out = v;
  return CRPX_OK;
}

void crpx_char_vector_free(crpx_char_vector_t vector)
{
  if (!vector) return;
  for (size_t i = 0; i < vector->size; ++i) free(vector->slots[i].text);
  free(vector->slots);
  free(vector);
}

size_t crpx_char_vector_size(const struct crpx_char_vector_struct *vector) { return vector ? vector->size : 0; }
const char *crpx_char_vector_get(const struct crpx_char_vector_struct *vector, size_t position)
{
  if (!vector || position >= vector->size) return NULL;
  return vector->slots[position].text ? vector->slots[position].text : "";
}
size_t crpx_char_vector_length(const struct crpx_char_vector_struct *vector, size_t position)
{
  return vector && position < vector->size ? vector->slots[position].length : 0;
}

crpx_status crpx_char_vector_resize(crpx_char_vector_t vector, size_t size)
{
  crpx_status status;
  if (!vector) return CRPX_ERROR_INVALID_ARGUMENT;
  status = reserve_slots(vector, size);
  if (status != CRPX_OK) return status;
  for (size_t i = size; i < vector->size; ++i) {
    free(vector->slots[i].text);
    memset(vector->slots + i, 0, sizeof(*vector->slots));
  }
  vector->size = size;
  return CRPX_OK;
}

crpx_status crpx_char_vector_set(crpx_char_vector_t vector, size_t position, const char *string)
{
  char *text;
  if (!vector || !string) return CRPX_ERROR_INVALID_ARGUMENT;
  if (position >= vector->size) return CRPX_ERROR_OUT_OF_RANGE;
  size_t length = strlen(string);
  if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  text = malloc(length + 1);
  if (!text) return CRPX_ERROR_ALLOCATION;
  memcpy(text, string, length + 1);
  free(vector->slots[position].text);
  vector->slots[position] = (struct string_slot){text, length, length + 1};
  return CRPX_OK;
}

crpx_status crpx_char_vector_append(crpx_char_vector_t vector, size_t position, const char *suffix)
{
  struct string_slot *slot;
  char *text;
  if (!vector || !suffix) return CRPX_ERROR_INVALID_ARGUMENT;
  if (position >= vector->size) return CRPX_ERROR_OUT_OF_RANGE;
  slot = vector->slots + position;
  size_t extra = strlen(suffix);
  if (!extra) return CRPX_OK;
  if (extra >= SIZE_MAX - slot->length) return CRPX_ERROR_OVERFLOW;
  size_t required = slot->length + extra + 1;
  /* Allocate separately on growth, so suffix may point inside any current slot. */
  if (required > slot->capacity) {
    size_t capacity = slot->capacity <= SIZE_MAX / 2 ? slot->capacity * 2 : required;
    if (capacity < required) capacity = required;
    text = malloc(capacity);
    if (!text) return CRPX_ERROR_ALLOCATION;
    if (slot->length) memcpy(text, slot->text, slot->length);
    memcpy(text + slot->length, suffix, extra + 1);
    free(slot->text);
    slot->text = text;
    slot->capacity = capacity;
  } else memmove(slot->text + slot->length, suffix, extra + 1);
  slot->length += extra;
  return CRPX_OK;
}

crpx_status crpx_char_vector_push(crpx_char_vector_t vector, const char *string)
{
  crpx_status status;
  char *copy;
  if (!vector || !string) return CRPX_ERROR_INVALID_ARGUMENT;
  if (vector->size == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  size_t length = strlen(string);
  if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  copy = malloc(length + 1);
  if (!copy) return CRPX_ERROR_ALLOCATION;
  memcpy(copy, string, length + 1);
  status = reserve_slots(vector, vector->size + 1);
  if (status != CRPX_OK) { free(copy); return status; }
  vector->slots[vector->size++] = (struct string_slot){copy, length, length + 1};
  return CRPX_OK;
}

crpx_status crpx_char_vector_clone(const struct crpx_char_vector_struct *source, crpx_char_vector_t *out)
{
  crpx_char_vector_t copy;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!source) return CRPX_ERROR_INVALID_ARGUMENT;
  status = crpx_char_vector_new(source->size, &copy);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < source->size; ++i) {
    status = crpx_char_vector_set(copy, i, crpx_char_vector_get(source, i));
    if (status != CRPX_OK) { crpx_char_vector_free(copy); return status; }
  }
  *out = copy;
  return CRPX_OK;
}

crpx_status crpx_char_vector_select(crpx_char_vector_t vector, const size_t *indices, size_t count)
{
  crpx_char_vector_t copy;
  crpx_status status;
  if (!vector || (!indices && count)) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < count; ++i) if (indices[i] >= vector->size) return CRPX_ERROR_OUT_OF_RANGE;
  status = crpx_char_vector_new(count, &copy);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < count; ++i) {
    status = crpx_char_vector_set(copy, i, crpx_char_vector_get(vector, indices[i]));
    if (status != CRPX_OK) { crpx_char_vector_free(copy); return status; }
  }
  struct crpx_char_vector_struct old = *vector;
  *vector = *copy;
  *copy = old;
  crpx_char_vector_free(copy);
  return CRPX_OK;
}

struct sort_slot { struct string_slot slot; size_t original; };
static int compare_lexical(const void *left, const void *right)
{
  const struct sort_slot *a = left, *b = right;
  int order = strcmp(a->slot.text ? a->slot.text : "", b->slot.text ? b->slot.text : "");
  return order ? order : (a->original > b->original) - (a->original < b->original);
}
static int compare_length(const void *left, const void *right)
{
  const struct sort_slot *a = left, *b = right;
  if (a->slot.length != b->slot.length) return a->slot.length < b->slot.length ? 1 : -1;
  return (a->original > b->original) - (a->original < b->original);
}

crpx_status crpx_char_vector_sort(crpx_char_vector_t vector, bool lexical, size_t *old_indices)
{
  struct sort_slot *order;
  if (!vector) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!vector->size) return CRPX_OK;
  if (vector->size > SIZE_MAX / sizeof(*order)) return CRPX_ERROR_OVERFLOW;
  order = malloc(vector->size * sizeof(*order));
  if (!order) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < vector->size; ++i) order[i] = (struct sort_slot){vector->slots[i], i};
  qsort(order, vector->size, sizeof(*order), lexical ? compare_lexical : compare_length);
  for (size_t i = 0; i < vector->size; ++i) {
    vector->slots[i] = order[i].slot;
    if (old_indices) old_indices[i] = order[i].original;
  }
  free(order);
  return CRPX_OK;
}

crpx_status crpx_char_vector_remove_empty(crpx_char_vector_t vector, size_t *removed)
{
  size_t count = 0;
  if (!vector) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < vector->size; ++i) {
    if (!vector->slots[i].length) free(vector->slots[i].text);
    else vector->slots[count++] = vector->slots[i];
  }
  if (removed) *removed = vector->size - count;
  if (vector->size > count) memset(vector->slots + count, 0, (vector->size - count) * sizeof(*vector->slots));
  vector->size = count;
  return CRPX_OK;
}

crpx_status crpx_char_vector_deduplicate(crpx_char_vector_t vector, size_t *removed)
{
  struct sort_slot *order;
  bool *duplicate;
  size_t count = 0;
  if (!vector) return CRPX_ERROR_INVALID_ARGUMENT;
  if (!vector->size) { if (removed) *removed = 0; return CRPX_OK; }
  if (vector->size > SIZE_MAX / sizeof(*order)) return CRPX_ERROR_OVERFLOW;
  order = malloc(vector->size * sizeof(*order));
  duplicate = calloc(vector->size, sizeof(*duplicate));
  if (!order || !duplicate) { free(order); free(duplicate); return CRPX_ERROR_ALLOCATION; }
  for (size_t i = 0; i < vector->size; ++i) order[i] = (struct sort_slot){vector->slots[i], i};
  qsort(order, vector->size, sizeof(*order), compare_lexical);
  for (size_t i = 1; i < vector->size; ++i) {
    const char *a = order[i-1].slot.text ? order[i-1].slot.text : "";
    const char *b = order[i].slot.text ? order[i].slot.text : "";
    if (!strcmp(a, b)) duplicate[order[i].original] = true;
  }
  for (size_t i = 0; i < vector->size; ++i) {
    if (duplicate[i]) free(vector->slots[i].text);
    else vector->slots[count++] = vector->slots[i];
  }
  if (removed) *removed = vector->size - count;
  if (vector->size > count) memset(vector->slots + count, 0, (vector->size - count) * sizeof(*vector->slots));
  vector->size = count;
  free(order);
  free(duplicate);
  return CRPX_OK;
}

bool crpx_char_vector_equal(const struct crpx_char_vector_struct *a, const struct crpx_char_vector_struct *b)
{
  if (!a || !b || a->size != b->size) return false;
  for (size_t i = 0; i < a->size; ++i)
    if (strcmp(crpx_char_vector_get(a, i), crpx_char_vector_get(b, i))) return false;
  return true;
}

crpx_status crpx_char_vector_find(const struct crpx_char_vector_struct *vector, const char *string, size_t *position)
{
  if (!vector || !string || !position) return CRPX_ERROR_INVALID_ARGUMENT;
  for (size_t i = 0; i < vector->size; ++i) if (!strcmp(crpx_char_vector_get(vector, i), string)) {
    *position = i;
    return CRPX_OK;
  }
  return CRPX_ERROR_NOT_FOUND;
}
