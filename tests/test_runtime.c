/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "global/global_variable.h"
#include "index_arrangement.h"
#include "quasi_random.h"
#include <pthread.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

enum { WORKERS = 4, DRAWS = 10000 };
struct worker {
  crpx_global_t context;
  uint64_t checksum;
  uint64_t identifier;
};

static void *shared_draws(void *argument)
{
  struct worker *worker = argument;
  worker->checksum = 0;
  for (int i = 0; i < DRAWS; i++) worker->checksum ^= crpx_random_64bits(worker->context);
  crpx_global_finalise(worker->context);
  return NULL;
}

static void *stream_draws(void *argument)
{
  struct worker *worker = argument;
  crpx_rng_stream stream;
  crpx_rng_stream_init(&stream, 0, worker->identifier);
  worker->checksum = 0;
  for (int i = 0; i < DRAWS; i++) worker->checksum ^= crpx_rng_stream_next(&stream);
  return NULL;
}

static int test_random(void)
{
  uint64_t star_state[4] = {1, 2, 3, 4}, plus_state[4] = {1, 2, 3, 4};
  CHECK(crpx_xoroshiro_star_seed256(star_state) == UINT64_C(11520));
  CHECK(crpx_xoroshiro_pp_seed256(plus_state) == UINT64_C(41943041));
  crpx_global_t a = crpx_global_init(123456789, NULL);
  crpx_global_t b = crpx_global_init(123456789, "fatal");
  CHECK(a && b);
  for (unsigned generator = 0; generator <= 20; generator++) {
    CHECK(crpx_set_random_generator_checked(a, generator, 123456789) == CRPX_OK);
    CHECK(crpx_set_random_generator_checked(b, generator, 123456789) == CRPX_OK);
    for (int i = 0; i < 500; i++) CHECK(crpx_random_64bits(a) == crpx_random_64bits(b));
  }
  CHECK(crpx_set_random_generator_checked(a, 255, 123) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(crpx_set_random_generator_checked(a, 0, 123) == CRPX_OK);
  CHECK(crpx_set_random_generator_checked(b, 0, 124) == CRPX_OK);
  CHECK(crpx_random_64bits(a) != crpx_random_64bits(b));
  uint64_t output = UINT64_MAX;
  CHECK(crpx_random_range_checked(a, 0, &output) == CRPX_ERROR_INVALID_ARGUMENT);
  CHECK(output == UINT64_MAX);
  CHECK(crpx_random_range_checked(a, 1, &output) == CRPX_OK && output == 0);
  for (int i = 0; i < 10000; i++) {
    CHECK(crpx_random_range_checked(a, UINT64_MAX, &output) == CRPX_OK && output < UINT64_MAX);
    CHECK(crpx_random_range_checked(a, 7, &output) == CRPX_OK && output < 7);
    double x = crpx_random_double(a);
    CHECK(x >= 0 && x < 1);
    CHECK(isfinite(crpx_random_normal(a, NULL)));
  }
  crpx_logger_set_level(a, CRPX_LOGLEVEL_FATAL);
  crpx_logger_error(a, "suppressed error still records failure");
  CHECK(a->error == 1);
  crpx_logger_message(255, __FILE__, __LINE__, a, "invalid level");
  CHECK(crpx_reallocarray(a, NULL, SIZE_MAX, 2) == NULL);
  crpx_global_finalise(a);
  crpx_global_finalise(b);
  return 0;
}

static int test_threads(void)
{
  crpx_global_t parallel = crpx_global_init(777, "fatal");
  crpx_global_t serial = crpx_global_init(777, "fatal");
  CHECK(parallel && serial);
  pthread_t threads[WORKERS];
  struct worker workers[WORKERS];
  for (uint64_t i = 0; i < WORKERS; i++) {
    workers[i].context = crpx_global_retain(parallel);
    workers[i].identifier = i;
    CHECK(pthread_create(&threads[i], NULL, shared_draws, &workers[i]) == 0);
  }
  uint64_t actual = 0, expected = 0;
  for (int i = 0; i < WORKERS; i++) {
    CHECK(pthread_join(threads[i], NULL) == 0);
    actual ^= workers[i].checksum;
  }
  for (int i = 0; i < WORKERS * DRAWS; i++) expected ^= crpx_random_64bits(serial);
  CHECK(actual == expected);
  for (int i = 0; i < WORKERS; i++) CHECK(pthread_create(&threads[i], NULL, stream_draws, &workers[i]) == 0);
  for (int i = 0; i < WORKERS; i++) {
    CHECK(pthread_join(threads[i], NULL) == 0);
    struct worker reference = { NULL, 0, (uint64_t)i };
    stream_draws(&reference);
    CHECK(reference.checksum == workers[i].checksum);
  }
  crpx_rng_stream a, b;
  crpx_rng_stream_init(&a, 0, 0);
  crpx_rng_stream_init(&b, 0, 1);
  CHECK(crpx_rng_stream_next(&a) != crpx_rng_stream_next(&b));
  uint64_t value = 99;
  CHECK(crpx_rng_stream_range(&a, 0, &value) == CRPX_ERROR_INVALID_ARGUMENT && value == 99);
  CHECK(crpx_rng_stream_range(&a, 1, &value) == CRPX_OK && value == 0);
  crpx_global_finalise(parallel);
  crpx_global_finalise(serial);
  return 0;
}

static int test_arrangements(void)
{
  crpx_global_t context = crpx_global_init(123, "fatal");
  CHECK(context);
  crpx_index_permutation_t permutation = crpx_index_permutation_new(context, 4);
  crpx_index_combination_t combination = crpx_index_combination_new(context, 5, 3);
  CHECK(permutation && combination);
  CHECK(crpx_index_combination_new(context, 2, 3) == NULL);
  CHECK(crpx_index_permutation_new(context, SIZE_MAX) == NULL);
  crpx_global_finalise(context);
  size_t count = 1;
  while (crpx_index_permutation_next(permutation)) { CHECK(crpx_index_permutation_is_valid(permutation)); count++; }
  CHECK(count == 24);
  while (crpx_index_permutation_prev(permutation)) count--;
  CHECK(count == 1);
  CHECK(crpx_index_permutation_swap(permutation, 0, 2));
  CHECK(crpx_index_permutation_inverse(permutation, permutation));
  CHECK(permutation->idx[0] == 2 && permutation->idx[2] == 0);
  CHECK(crpx_index_permutation_combine(permutation, permutation, permutation));
  for (size_t i = 0; i < 4; i++) CHECK(permutation->idx[i] == i);
  count = 1;
  while (crpx_index_combination_next(combination)) { CHECK(crpx_index_combination_is_valid(combination)); count++; }
  CHECK(count == 10);
  while (crpx_index_combination_prev(combination)) count--;
  CHECK(count == 1);
  del_crpx_index_combination(combination);
  del_crpx_index_permutation(permutation);
  crpx_index_permutation_t empty = crpx_index_permutation_new(NULL, 0);
  crpx_index_combination_t empty_combination = crpx_index_combination_new(NULL, 0, 0);
  CHECK(empty && empty_combination);
  CHECK(crpx_index_permutation_is_valid(empty));
  CHECK(crpx_index_combination_is_valid(empty_combination));
  del_crpx_index_permutation(empty);
  del_crpx_index_combination(empty_combination);
  return 0;
}

static int test_quasi_and_bits(void)
{
  crpx_global_t context = crpx_global_init(321, "fatal");
  CHECK(context);
  CHECK(new_crpx_quasi_random(context, 0) == NULL);
  for (size_t dimension = 1; dimension <= 1230; dimension += dimension < 3 ? 1 : 1227) {
    crpx_quasi_random_t q = new_crpx_quasi_random(context, dimension);
    CHECK(q);
    for (int i = 0; i < 20; i++) {
      crpx_quasi_random_next_halton(q);
      for (size_t j = 0; j < dimension; j++) CHECK(q->r[j] >= 0 && q->r[j] < 1);
      crpx_quasi_random_next_halton_original(q);
      for (size_t j = 0; j < dimension; j++) CHECK(q->r[j] >= 0 && q->r[j] < 1);
    }
    del_crpx_quasi_random(q);
  }
  crpx_quasi_random_t retained = new_crpx_quasi_random(context, 3);
  CHECK(retained);
  crpx_global_finalise(context);
  crpx_quasi_random_next_korobov(retained);
  del_crpx_quasi_random(retained);
  const uint64_t values[] = {0, 1, UINT64_MAX, UINT64_C(0x12345678abcdef90)};
  for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++)
    CHECK(crpx_deinterleave_64bits(crpx_interleave_64bits(values[i])) == values[i]);
  CHECK(crpx_choose_n_k(7, 3) == 35);
  CHECK(crpx_choose_n_k(100000, 2) == 0);
  CHECK(crpx_choose_n_k(-1, 1) == 0);
  int combination[3] = {0};
  crpx_ordered_combination_n_k(combination, 5, 1, 4);
  CHECK(combination[0] == 5);
  crpx_ordered_combination_n_k(combination, 5, 3, 9);
  CHECK(combination[0] == 3 && combination[1] == 4 && combination[2] == 5);
  return 0;
}

int main(void)
{
  CHECK(test_random() == 0);
  CHECK(test_threads() == 0);
  CHECK(test_arrangements() == 0);
  CHECK(test_quasi_and_bits() == 0);
  return 0;
}
