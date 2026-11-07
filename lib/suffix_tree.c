/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "suffix_tree.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct crpx_suffix_index {
  unsigned char *data;
  size_t length, *suffixes, *lcp;
};
struct suffix_record { size_t position, first, second; };

static int compare_records(const void *left, const void *right)
{
  const struct suffix_record *a = left, *b = right;
  if (a->first != b->first) return a->first < b->first ? -1 : 1;
  if (a->second != b->second) return a->second < b->second ? -1 : 1;
  return 0;
}

void crpx_suffix_index_free(crpx_suffix_index *index)
{
  if (!index) return;
  free(index->data); free(index->suffixes); free(index->lcp); free(index);
}

size_t crpx_suffix_index_length(const crpx_suffix_index *index) { return index ? index->length : 0; }

crpx_status crpx_suffix_index_new(const void *data, size_t length, crpx_suffix_index **result)
{
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!data && length) return CRPX_ERROR_INVALID_ARGUMENT;
  if (length > SIZE_MAX / sizeof(struct suffix_record)) return CRPX_ERROR_OVERFLOW;
  crpx_suffix_index *index = calloc(1, sizeof(*index));
  if (!index) return CRPX_ERROR_ALLOCATION;
  index->length = length;
  if (!length) { *result = index; return CRPX_OK; }
  index->data = malloc(length);
  index->suffixes = malloc(length * sizeof(size_t));
  index->lcp = calloc(length, sizeof(size_t));
  size_t *ranks = malloc(length * sizeof(size_t));
  size_t *next_ranks = malloc(length * sizeof(size_t));
  struct suffix_record *records = malloc(length * sizeof(*records));
  if (!index->data || !index->suffixes || !index->lcp || !ranks || !next_ranks || !records) {
    free(ranks); free(next_ranks); free(records); crpx_suffix_index_free(index);
    return CRPX_ERROR_ALLOCATION;
  }
  memcpy(index->data, data, length);
  for (size_t i = 0; i < length; i++) { ranks[i] = (size_t)index->data[i] + 1; index->suffixes[i] = i; }
  for (size_t width = 1; width < length;) {
    for (size_t i = 0; i < length; i++) {
      records[i].position = i;
      records[i].first = ranks[i];
      records[i].second = width < length - i ? ranks[i + width] : 0;
    }
    qsort(records, length, sizeof(*records), compare_records);
    size_t rank = 1;
    for (size_t i = 0; i < length; i++) {
      if (i && compare_records(&records[i - 1], &records[i])) rank++;
      next_ranks[records[i].position] = rank;
      index->suffixes[i] = records[i].position;
    }
    size_t *temporary = ranks; ranks = next_ranks; next_ranks = temporary;
    if (rank == length) break;
    width = width > length / 2 ? length : width * 2;
  }
  /* Kasai's linear-time adjacent LCP construction. */
  for (size_t i = 0; i < length; i++) ranks[index->suffixes[i]] = i;
  size_t common = 0;
  for (size_t i = 0; i < length; i++) {
    size_t rank = ranks[i];
    if (!rank) { common = 0; continue; }
    size_t previous = index->suffixes[rank - 1];
    while (common < length - i && common < length - previous && index->data[i + common] == index->data[previous + common])
      common++;
    index->lcp[rank] = common;
    if (common) common--;
  }
  free(ranks); free(next_ranks); free(records);
  *result = index;
  return CRPX_OK;
}

static int compare_pattern(const crpx_suffix_index *index, size_t rank, const unsigned char *pattern, size_t length)
{
  size_t position = index->suffixes[rank];
  size_t available = index->length - position;
  size_t common = length < available ? length : available;
  int comparison = memcmp(index->data + position, pattern, common);
  if (comparison) return comparison;
  return available < length ? -1 : 0;
}

crpx_status crpx_suffix_index_find(const crpx_suffix_index *index, const void *pattern, size_t length,
                                 size_t *first, size_t *count)
{
  if (!index || !pattern || !length || !first || !count) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t lower = 0, upper = index->length;
  while (lower < upper) {
    size_t middle = lower + (upper - lower) / 2;
    if (compare_pattern(index, middle, pattern, length) < 0) lower = middle + 1;
    else upper = middle;
  }
  size_t start = lower;
  upper = index->length;
  while (lower < upper) {
    size_t middle = lower + (upper - lower) / 2;
    if (compare_pattern(index, middle, pattern, length) <= 0) lower = middle + 1;
    else upper = middle;
  }
  if (lower == start) return CRPX_ERROR_NOT_FOUND;
  *first = start; *count = lower - start;
  return CRPX_OK;
}

crpx_status crpx_suffix_index_position(const crpx_suffix_index *index, size_t suffix_rank, size_t *position)
{
  if (!index || !position) return CRPX_ERROR_INVALID_ARGUMENT;
  if (suffix_rank >= index->length) return CRPX_ERROR_OUT_OF_RANGE;
  *position = index->suffixes[suffix_rank];
  return CRPX_OK;
}

crpx_status crpx_suffix_index_lcp(const crpx_suffix_index *index, size_t first, size_t second, size_t *length)
{
  if (!index || !length) return CRPX_ERROR_INVALID_ARGUMENT;
  if (first > index->length || second > index->length) return CRPX_ERROR_OUT_OF_RANGE;
  size_t common = 0;
  while (common < index->length - first && common < index->length - second &&
         index->data[first + common] == index->data[second + common]) common++;
  *length = common;
  return CRPX_OK;
}

crpx_status crpx_suffix_index_longest_repeat(const crpx_suffix_index *index, size_t *position, size_t *length)
{
  if (!index || !position || !length) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t longest = 0, coordinate = 0;
  for (size_t i = 1; i < index->length; i++) {
    if (index->lcp[i] > longest) { longest = index->lcp[i]; coordinate = index->suffixes[i]; }
  }
  *position = coordinate; *length = longest;
  return CRPX_OK;
}
