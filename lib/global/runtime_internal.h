/* Internal synchronisation, deliberately absent from the installed API. */
#ifndef CRPX_RUNTIME_INTERNAL_H
#define CRPX_RUNTIME_INTERNAL_H
#include <stdatomic.h>
struct crpx_runtime_state {
  atomic_flag rng_lock;
  atomic_flag logger_lock;
  atomic_uint references;
};
static inline void crpx_runtime_lock (atomic_flag *lock)
{
  while (atomic_flag_test_and_set_explicit(lock, memory_order_acquire)) { }
}
static inline void crpx_runtime_unlock (atomic_flag *lock)
{
  atomic_flag_clear_explicit(lock, memory_order_release);
}
#endif
