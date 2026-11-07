/* Reimplementation of biomcmc-lib hash containers (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "hashtable.h"
#include <stdlib.h>
#include <string.h>

struct hash_entry { unsigned char *key; size_t length; uint64_t hash, value; unsigned char state; };
struct crpx_hashtable_struct { struct hash_entry *entries; size_t capacity, size, occupied; uint64_t seed; };

static uint64_t hash_key(const void *key, size_t length, uint64_t seed)
{
  const unsigned char *bytes = key;
  uint64_t hash = UINT64_C(14695981039346656037) ^ seed;
  for (size_t i = 0; i < length; ++i) hash = (hash ^ bytes[i]) * UINT64_C(1099511628211);
  hash ^= hash >> 30;
  hash *= UINT64_C(0xbf58476d1ce4e5b9);
  hash ^= hash >> 27;
  hash *= UINT64_C(0x94d049bb133111eb);
  return hash ^ (hash >> 31);
}

static size_t locate(const struct crpx_hashtable_struct *table, const void *key, size_t length, uint64_t hash,
                     bool *found)
{
  size_t position = (size_t) hash & (table->capacity - 1), first_deleted = SIZE_MAX;
  for (size_t i = 0; i < table->capacity; ++i, position = (position + 1) & (table->capacity - 1)) {
    const struct hash_entry *entry = table->entries + position;
    if (!entry->state) { *found = false; return first_deleted == SIZE_MAX ? position : first_deleted; }
    if (entry->state == 2) { if (first_deleted == SIZE_MAX) first_deleted = position; continue; }
    if (entry->hash == hash && entry->length == length && (!length || !memcmp(entry->key, key, length))) {
      *found = true;
      return position;
    }
  }
  *found = false;
  return first_deleted;
}

static crpx_status rehash(crpx_hashtable_t table, size_t capacity)
{
  struct hash_entry *entries;
  struct crpx_hashtable_struct fresh = {0};
  if (capacity > SIZE_MAX / sizeof(*entries)) return CRPX_ERROR_OVERFLOW;
  entries = calloc(capacity, sizeof(*entries));
  if (!entries) return CRPX_ERROR_ALLOCATION;
  fresh.entries = entries;
  fresh.capacity = capacity;
  for (size_t i = 0; i < table->capacity; ++i) if (table->entries[i].state == 1) {
    bool found;
    size_t position = locate(&fresh, table->entries[i].key, table->entries[i].length, table->entries[i].hash, &found);
    entries[position] = table->entries[i];
  }
  free(table->entries);
  table->entries = entries;
  table->capacity = capacity;
  table->occupied = table->size;
  return CRPX_OK;
}

crpx_status crpx_hashtable_new(size_t expected_items, uint64_t seed, crpx_hashtable_t *out)
{
  size_t capacity = 16;
  crpx_hashtable_t table;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (expected_items > SIZE_MAX / 2) return CRPX_ERROR_OVERFLOW;
  while (capacity / 2 < expected_items) {
    if (capacity > SIZE_MAX / 2) return CRPX_ERROR_OVERFLOW;
    capacity *= 2;
  }
  table = calloc(1, sizeof(*table));
  if (!table) return CRPX_ERROR_ALLOCATION;
  status = rehash(table, capacity);
  if (status != CRPX_OK) { free(table); return status; }
  table->seed = seed;
  *out = table;
  return CRPX_OK;
}

void crpx_hashtable_clear(crpx_hashtable_t table)
{
  if (!table) return;
  for (size_t i = 0; i < table->capacity; ++i) if (table->entries[i].state == 1) free(table->entries[i].key);
  memset(table->entries, 0, table->capacity * sizeof(*table->entries));
  table->size = table->occupied = 0;
}

void crpx_hashtable_free(crpx_hashtable_t table)
{
  if (table) { crpx_hashtable_clear(table); free(table->entries); free(table); }
}

size_t crpx_hashtable_size(const struct crpx_hashtable_struct *table) { return table ? table->size : 0; }

crpx_status crpx_hashtable_put(crpx_hashtable_t table, const void *key, size_t length, uint64_t value, bool *replaced)
{
  crpx_status status;
  bool found;
  unsigned char *copy;
  if (!table || (!key && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t hash = hash_key(key, length, table->seed);
  size_t position = locate(table, key, length, hash, &found);
  if (found) {
    table->entries[position].value = value;
    if (replaced) *replaced = true;
    return CRPX_OK;
  }
  /* Copy before rehash so a key obtained by iteration can be inserted safely. */
  copy = malloc(length ? length : 1);
  if (!copy) return CRPX_ERROR_ALLOCATION;
  if (length) memcpy(copy, key, length);
  if (table->occupied >= table->capacity - table->capacity / 4) {
    size_t capacity = table->capacity;
    if (table->size >= table->capacity / 2) {
      if (capacity > SIZE_MAX / 2) { free(copy); return CRPX_ERROR_OVERFLOW; }
      capacity *= 2;
    }
    status = rehash(table, capacity);
    if (status != CRPX_OK) { free(copy); return status; }
    position = locate(table, copy, length, hash, &found);
  }
  if (position == SIZE_MAX) { free(copy); return CRPX_ERROR_OVERFLOW; }
  if (!table->entries[position].state) ++table->occupied;
  table->entries[position] = (struct hash_entry){copy, length, hash, value, 1};
  ++table->size;
  if (replaced) *replaced = false;
  return CRPX_OK;
}

crpx_status crpx_hashtable_get(const struct crpx_hashtable_struct *table, const void *key, size_t length, uint64_t *value)
{
  bool found;
  if (!table || (!key && length) || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t position = locate(table, key, length, hash_key(key, length, table->seed), &found);
  if (!found) return CRPX_ERROR_NOT_FOUND;
  *value = table->entries[position].value;
  return CRPX_OK;
}

crpx_status crpx_hashtable_remove(crpx_hashtable_t table, const void *key, size_t length)
{
  bool found;
  if (!table || (!key && length)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t position = locate(table, key, length, hash_key(key, length, table->seed), &found);
  if (!found) return CRPX_ERROR_NOT_FOUND;
  free(table->entries[position].key);
  table->entries[position].key = NULL;
  table->entries[position].state = 2;
  --table->size;
  return CRPX_OK;
}

crpx_status crpx_hashtable_next(const struct crpx_hashtable_struct *table, size_t *cursor,
                              const void **key, size_t *length, uint64_t *value)
{
  if (!table || !cursor || !key || !length || !value) return CRPX_ERROR_INVALID_ARGUMENT;
  while (*cursor < table->capacity) {
    const struct hash_entry *entry = table->entries + (*cursor)++;
    if (entry->state == 1) { *key = entry->key; *length = entry->length; *value = entry->value; return CRPX_OK; }
  }
  return CRPX_ERROR_NOT_FOUND;
}

crpx_status crpx_hashtable_put_string(crpx_hashtable_t table, const char *key, uint64_t value, bool *replaced)
{
  return key ? crpx_hashtable_put(table, key, strlen(key), value, replaced) : CRPX_ERROR_INVALID_ARGUMENT;
}
crpx_status crpx_hashtable_get_string(const struct crpx_hashtable_struct *table, const char *key, uint64_t *value)
{
  return key ? crpx_hashtable_get(table, key, strlen(key), value) : CRPX_ERROR_INVALID_ARGUMENT;
}
