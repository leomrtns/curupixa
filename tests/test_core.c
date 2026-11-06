/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "bipartition.h"
#include "char_vector.h"
#include "distance_matrix.h"
#include "hashtable.h"

#define OK(expression) assert((expression) == CRPX_OK)

static void test_bits(void)
{
  const size_t sizes[] = {0, 1, 63, 64, 65, 127, 128, 129};
  for (size_t n = 0; n < sizeof(sizes) / sizeof(*sizes); ++n) {
    crpx_bipartition_t a, b, c;
    size_t size = sizes[n], count = 0;
    OK(crpx_bipartition_new(size, &a));
    OK(crpx_bipartition_new(size, &b));
    OK(crpx_bipartition_new(size, &c));
    for (size_t i = 0; i < size; ++i) if (i % 3 == 0) { OK(crpx_bipartition_set(a, i, true)); ++count; }
    assert(crpx_bipartition_count(a) == count);
    OK(crpx_bipartition_complement(b, a));
    assert(crpx_bipartition_count(b) == size - count);
    assert(crpx_bipartition_same_split(a, b));
    OK(crpx_bipartition_combine(c, a, b, CRPX_BITS_OR));
    assert(crpx_bipartition_count(c) == size);
    OK(crpx_bipartition_combine(c, a, b, CRPX_BITS_AND));
    assert(crpx_bipartition_count(c) == 0);
    OK(crpx_bipartition_combine(b, b, a, CRPX_BITS_XORNOT));
    assert(crpx_bipartition_count(b) == 0);
    OK(crpx_bipartition_resize(a, size + 70));
    assert(crpx_bipartition_count(a) == count);
    OK(crpx_bipartition_resize(a, size));
    assert(crpx_bipartition_count(a) == count);
    crpx_bipartition_free(a); crpx_bipartition_free(b); crpx_bipartition_free(c);
  }
  crpx_bipartition_t a, b;
  OK(crpx_bipartition_new(4, &a));
  OK(crpx_bipartition_new(4, &b));
  OK(crpx_bipartition_set(a, 0, true)); OK(crpx_bipartition_set(a, 1, true));
  OK(crpx_bipartition_set(b, 0, true)); OK(crpx_bipartition_set(b, 2, true));
  assert(!crpx_bipartition_compatible(a, b));
  crpx_bipartition_canonicalise(a);
  bool value;
  OK(crpx_bipartition_get(a, 0, &value)); assert(!value);
  size_t indices[2], count;
  assert(crpx_bipartition_indices(a, NULL, 0, &count) == CRPX_ERROR_OUT_OF_RANGE && count == 2);
  OK(crpx_bipartition_indices(a, indices, 2, &count));
  assert(indices[0] == 2 && indices[1] == 3);
  assert(crpx_bipartition_set(a, 4, true) == CRPX_ERROR_OUT_OF_RANGE);
  crpx_bipartition_free(a); crpx_bipartition_free(b);
}

static void test_strings(void)
{
  crpx_char_vector_t v, copy;
  size_t removed, old[5];
  OK(crpx_char_vector_new(1, &v));
  assert(!strcmp(crpx_char_vector_get(v, 0), ""));
  OK(crpx_char_vector_set(v, 0, "taxon"));
  OK(crpx_char_vector_append(v, 0, crpx_char_vector_get(v, 0) + 2));
  assert(!strcmp(crpx_char_vector_get(v, 0), "taxonxon"));
  OK(crpx_char_vector_set(v, 0, crpx_char_vector_get(v, 0) + 5));
  OK(crpx_char_vector_push(v, "zebra")); OK(crpx_char_vector_push(v, ""));
  OK(crpx_char_vector_push(v, "xon")); OK(crpx_char_vector_push(v, "ant"));
  OK(crpx_char_vector_clone(v, &copy)); assert(crpx_char_vector_equal(v, copy));
  OK(crpx_char_vector_sort(v, true, old));
  assert(!strcmp(crpx_char_vector_get(v, 0), "") && old[0] == 2);
  OK(crpx_char_vector_remove_empty(v, &removed)); assert(removed == 1);
  OK(crpx_char_vector_deduplicate(v, &removed)); assert(removed == 1);
  assert(crpx_char_vector_size(v) == 3);
  const size_t selection[] = {2, 0, 2};
  OK(crpx_char_vector_select(v, selection, 3));
  assert(!strcmp(crpx_char_vector_get(v, 0), "zebra"));
  assert(!strcmp(crpx_char_vector_get(v, 1), "ant"));
  assert(crpx_char_vector_set(v, 10, "x") == CRPX_ERROR_OUT_OF_RANGE);
  OK(crpx_char_vector_resize(v, 0)); OK(crpx_char_vector_resize(v, 3));
  assert(!strcmp(crpx_char_vector_get(v, 2), ""));
  crpx_char_vector_free(v); crpx_char_vector_free(copy);
}

static void test_map(void)
{
  crpx_hashtable_t table;
  char key[32];
  uint64_t value;
  bool replaced;
  OK(crpx_hashtable_new(0, 1234, &table));
  for (size_t i = 0; i < 2000; ++i) {
    snprintf(key, sizeof(key), "taxon-%zu", i);
    OK(crpx_hashtable_put_string(table, key, i, &replaced)); assert(!replaced);
  }
  for (size_t i = 0; i < 2000; i += 2) {
    snprintf(key, sizeof(key), "taxon-%zu", i); OK(crpx_hashtable_remove(table, key, strlen(key)));
  }
  for (size_t i = 0; i < 2000; ++i) {
    snprintf(key, sizeof(key), "taxon-%zu", i);
    crpx_status status = crpx_hashtable_get_string(table, key, &value);
    assert(i % 2 ? status == CRPX_OK && value == i : status == CRPX_ERROR_NOT_FOUND);
  }
  const unsigned char binary[] = {'a', 0, 'b'};
  OK(crpx_hashtable_put(table, binary, 3, 42, NULL));
  OK(crpx_hashtable_put(table, binary, 3, 43, &replaced)); assert(replaced);
  OK(crpx_hashtable_get(table, binary, 3, &value)); assert(value == 43);
  OK(crpx_hashtable_put(table, NULL, 0, 7, NULL));
  OK(crpx_hashtable_get(table, NULL, 0, &value)); assert(value == 7);
  size_t cursor = 0, length, seen = 0;
  const void *bytes;
  while (crpx_hashtable_next(table, &cursor, &bytes, &length, &value) == CRPX_OK) ++seen;
  assert(seen == crpx_hashtable_size(table));
  crpx_hashtable_clear(table); assert(!crpx_hashtable_size(table));
  crpx_hashtable_free(table);
}

static void test_distances(void)
{
  crpx_distance_matrix_t genes, output;
  crpx_species_distances_t species;
  double value, minimum;
  uint64_t count;
  OK(crpx_distance_matrix_new(4, &genes));
  OK(crpx_distance_matrix_get(genes, 0, 1, &value)); assert(isnan(value));
  OK(crpx_distance_matrix_fill(genes, 2));
  OK(crpx_distance_matrix_set(genes, 0, 2, 4));
  OK(crpx_distance_matrix_get(genes, 2, 0, &value)); assert(value == 4);
  assert(crpx_distance_matrix_set(genes, 0, 1, -1) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(crpx_distance_matrix_new(SIZE_MAX, &output) == CRPX_ERROR_OVERFLOW && output == NULL);
  OK(crpx_species_distances_new(3, &species));
  const size_t map[] = {0, 0, 1, SIZE_MAX};
  OK(crpx_species_distances_add_genes(species, genes, map));
  OK(crpx_species_distances_get(species, 0, 1, &value, &minimum, &count));
  assert(value == 3 && minimum == 2 && count == 2);
  OK(crpx_species_distances_merge(species, species));
  OK(crpx_species_distances_get(species, 0, 1, &value, &minimum, &count)); assert(value == 3 && count == 4);
  OK(crpx_distance_matrix_new(3, &output));
  OK(crpx_species_distances_export(species, output, 1));
  OK(crpx_distance_matrix_get(output, 0, 1, &value)); assert(value == 2);
  OK(crpx_distance_matrix_get(output, 0, 2, &value)); assert(isnan(value));
  crpx_distance_matrix_free(genes); crpx_distance_matrix_free(output); crpx_species_distances_free(species);
}

int main(void)
{
  test_bits(); test_strings(); test_map(); test_distances();
  puts("core containers: passed");
  return 0;
}
