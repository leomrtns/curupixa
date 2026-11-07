/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "minimiser_simplex.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

static crpx_status rosenbrock(const double *x, size_t n, void *params, double *value)
{
  (void)params; assert(n == 2);
  *value = (1-x[0])*(1-x[0])+100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
  return CRPX_OK;
}
static crpx_status fail(const double *x, size_t n, void *params, double *value)
{
  (void)x; (void)n; (void)params; (void)value;
  return CRPX_ERROR_IO;
}
int main(void)
{
  crpx_simplex *s;
  const double point[] = {-1.2, 1}, steps[] = {0.5, 0.5};
  size_t iterations;
  assert(crpx_simplex_new(2, rosenbrock, NULL, &s) == CRPX_OK);
  assert(crpx_simplex_initialise(s, point, steps) == CRPX_OK);
  assert(crpx_simplex_run(s, 1, 1e-10, 1e-12, &iterations) == CRPX_ERROR_CONVERGENCE && iterations == 1);
  assert(crpx_simplex_run(s, 2000, 1e-10, 1e-12, &iterations) == CRPX_OK);
  const double *best = crpx_simplex_best_point(s);
  assert(fabs(best[0]-1) < 1e-7 && fabs(best[1]-1) < 1e-7 && crpx_simplex_best_value(s) < 1e-18);
  crpx_simplex_free(s);
  assert(crpx_simplex_new(2, fail, NULL, &s) == CRPX_OK);
  assert(crpx_simplex_initialise(s, point, steps) == CRPX_ERROR_IO);
  assert(!crpx_simplex_best_point(s));
  crpx_simplex_free(s);
  puts("bounded Nelder-Mead minimiser: passed");
  return 0;
}
