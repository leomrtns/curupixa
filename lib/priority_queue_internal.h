/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file priority_queue_internal.h
    Indexed minimum heap for graph searches and OPTICS. Equal priorities use the smallest item identifier. */
#ifndef CRPX_PRIORITY_QUEUE_INTERNAL_H
#define CRPX_PRIORITY_QUEUE_INTERNAL_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "status.h"

typedef struct {
  size_t size, capacity;
  size_t *items, *positions;
  long double *keys;
} crpx_priority_queue;

static inline void
crpx_queue_free (crpx_priority_queue *queue)
{
  free (queue->items); free (queue->positions); free (queue->keys);
  *queue = (crpx_priority_queue){0};
}

static inline crpx_status
crpx_queue_init (crpx_priority_queue *queue, size_t capacity)
{
  *queue = (crpx_priority_queue){0};
  if (capacity > SIZE_MAX / sizeof (long double) || capacity > SIZE_MAX / sizeof (size_t))
    return CRPX_ERROR_OVERFLOW;
  if (!capacity) return CRPX_OK;
  queue->items = malloc (capacity * sizeof (*queue->items));
  queue->positions = malloc (capacity * sizeof (*queue->positions));
  queue->keys = malloc (capacity * sizeof (*queue->keys));
  if (!queue->items || !queue->positions || !queue->keys) { crpx_queue_free (queue); return CRPX_ERROR_ALLOCATION; }
  queue->capacity = capacity;
  for (size_t i = 0; i < capacity; i++) queue->positions[i] = SIZE_MAX;
  return CRPX_OK;
}

static inline int
crpx_queue_before (const crpx_priority_queue *queue, size_t first, size_t second)
{
  return queue->keys[first] < queue->keys[second] ||
         (queue->keys[first] == queue->keys[second] && first < second);
}

static inline void
crpx_queue_offer (crpx_priority_queue *queue, size_t item, long double key)
{ // Insert an absent item or decrease its priority. The caller supplies valid IDs and ordered keys (finite or +INFINITY).
  size_t position = queue->positions[item];
  if (position == SIZE_MAX) position = queue->size++;
  else if (key >= queue->keys[item]) return;
  queue->keys[item] = key;
  while (position) {
    size_t parent = (position - 1) / 2, other = queue->items[parent];
    if (!crpx_queue_before (queue, item, other)) break;
    queue->items[position] = other; queue->positions[other] = position;
    position = parent;
  }
  queue->items[position] = item; queue->positions[item] = position;
}

static inline size_t
crpx_queue_pop (crpx_priority_queue *queue)
{
  if (!queue->size) return SIZE_MAX;
  size_t result = queue->items[0], item = queue->items[--queue->size], position = 0;
  queue->positions[result] = SIZE_MAX;
  if (!queue->size) return result;
  while (position < queue->size / 2) {
    size_t child = 2 * position + 1;
    if (child + 1 < queue->size && crpx_queue_before (queue, queue->items[child + 1], queue->items[child])) child++;
    if (!crpx_queue_before (queue, queue->items[child], item)) break;
    queue->items[position] = queue->items[child]; queue->positions[queue->items[position]] = position;
    position = child;
  }
  queue->items[position] = item; queue->positions[item] = position;
  return result;
}
#endif
