/* SPDX-License-Identifier: GPL-3.0-or-later
 * Independent generator/uniformisation/resolvent checks for the analytic models. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "substitution_model.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK(call) assert((call) == CRPX_OK)

static void close_value(double actual, long double expected, long double tolerance)
{
  long double error = fabsl((long double)actual-expected);
  if (error > tolerance*(fabsl(expected)+1e-290L)) {
    fprintf(stderr, "actual %.17g expected %.21Lg error %.6Lg tolerance %.6Lg\n",
            actual, expected, error, tolerance*(fabsl(expected)+1e-290L));
    assert(0);
  }
}

/* Build the generator directly from exchangeabilities and normalise its
 * stationary exit rate. No use of the implementation's returned generator. */
static void reference_generator(size_t n, const double *weights, double kappa, int hky,
                                long double *pi, long double *q)
{
  long double maximum = 0, sum = 0, mean = 0;
  for (size_t i = 0; i < n; ++i) if (weights[i] > maximum) maximum = weights[i];
  for (size_t i = 0; i < n; ++i) sum += (long double)weights[i]/maximum;
  for (size_t i = 0; i < n; ++i) pi[i] = ((long double)weights[i]/maximum)/sum;
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    int transition = hky && ((i == 0 && j == 2) || (i == 2 && j == 0) ||
                             (i == 1 && j == 3) || (i == 3 && j == 1));
    q[i*n+j] = i == j ? 0 : pi[j]*(transition ? kappa : 1);
    mean += pi[i]*q[i*n+j];
  }
  if (n == 1) { q[0] = 0; return; }
  for (size_t i = 0; i < n; ++i) {
    long double row = 0;
    for (size_t j = 0; j < n; ++j) if (i != j) { q[i*n+j] /= mean; row += q[i*n+j]; }
    q[i*n+i] = -row;
  }
}

static void multiply(const long double *a, const long double *b, long double *out, size_t n)
{
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    long double sum = 0;
    for (size_t k = 0; k < n; ++k) sum += a[i*n+k]*b[k*n+j];
    out[i*n+j] = sum;
  }
}

/* Uniformisation is a positive Poisson-weighted sum of powers of I+Q/nu.
 * It is independent of the analytic F81/HKY projectors under test. */
static void reference_exponential(const long double *q, size_t n, long double t, long double *out)
{
  long double nu = 0;
  for (size_t i = 0; i < n; ++i) if (-q[i*n+i] > nu) nu = -q[i*n+i];
  if (!nu || !t) {
    for (size_t i = 0; i < n*n; ++i) out[i] = i/n == i%n;
    return;
  }
  long double *r = malloc(n*n*sizeof(*r)), *power = calloc(n*n,sizeof(*power));
  long double *next = malloc(n*n*sizeof(*next));
  assert(r && power && next);
  for (size_t i = 0; i < n*n; ++i) {
    r[i] = q[i]/nu+(i/n == i%n);
    power[i] = i/n == i%n;
    out[i] = 0;
  }
  long double weight = expl(-nu*t);
  for (size_t step = 0;; ++step) {
    assert(step < 10000);
    for (size_t i = 0; i < n*n; ++i) out[i] += weight*power[i];
    if (step > nu*t+10 && weight < 1e-28L) break;
    multiply(power,r,next,n);
    memcpy(power,next,n*n*sizeof(*power));
    weight *= nu*t/(step+1);
  }
  free(r); free(power); free(next);
}

/* Solve (I-mu Q) X=I with pivoted Gauss-Jordan elimination. */
static void reference_resolvent(const long double *q, size_t n, long double mean, long double *out)
{
  long double *a = malloc(n*2*n*sizeof(*a)); assert(a);
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) {
    a[i*2*n+j] = (i == j)-mean*q[i*n+j];
    a[i*2*n+n+j] = i == j;
  }
  for (size_t k = 0; k < n; ++k) {
    size_t pivot = k;
    for (size_t i = k+1; i < n; ++i) if (fabsl(a[i*2*n+k]) > fabsl(a[pivot*2*n+k])) pivot = i;
    assert(a[pivot*2*n+k] != 0);
    for (size_t j = 0; j < 2*n; ++j) {
      long double tmp = a[k*2*n+j]; a[k*2*n+j] = a[pivot*2*n+j]; a[pivot*2*n+j] = tmp;
    }
    long double divisor = a[k*2*n+k];
    for (size_t j = 0; j < 2*n; ++j) a[k*2*n+j] /= divisor;
    for (size_t i = 0; i < n; ++i) if (i != k) {
      long double factor = a[i*2*n+k];
      for (size_t j = 0; j < 2*n; ++j) a[i*2*n+j] -= factor*a[k*2*n+j];
    }
  }
  for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) out[i*n+j] = a[i*2*n+n+j];
  free(a);
}

static void matrix_invariants(const double *p, const long double *pi, size_t n)
{
  for (size_t i = 0; i < n; ++i) {
    long double sum = 0;
    for (size_t j = 0; j < n; ++j) {
      assert(isfinite(p[i*n+j]) && p[i*n+j] >= 0 && p[i*n+j] <= 1+1e-14);
      sum += p[i*n+j];
      long double lhs = pi[i]*p[i*n+j], rhs = pi[j]*p[j*n+i];
      assert(fabsl(lhs-rhs) <= 2e-13L*(fabsl(lhs)+fabsl(rhs)+1e-300L));
    }
    assert(fabsl(sum-1) < 2e-13L);
  }
}

static void compare_model(size_t n, const double *frequencies, double kappa, int hky)
{
  crpx_substitution_model *m;
  if (hky) OK(crpx_substitution_model_hky(frequencies,kappa,&m));
  else OK(crpx_substitution_model_f81(n,frequencies,&m));
  long double *q = malloc(n*n*sizeof(*q)), *pi = malloc(n*sizeof(*pi)), *reference = malloc(n*n*sizeof(*reference));
  double *actual = malloc(n*n*sizeof(*actual)), *generator = malloc(n*n*sizeof(*generator));
  double *child = malloc(n*sizeof(*child)), *applied = malloc(n*sizeof(*applied));
  assert(q && pi && reference && actual && generator && child && applied);
  reference_generator(n,frequencies,kappa,hky,pi,q);
  OK(crpx_substitution_model_generator(m,generator,n*n));
  long double mean = 0;
  for (size_t i = 0; i < n*n; ++i) close_value(generator[i],q[i],3e-13L);
  for (size_t i = 0; i < n; ++i) { mean -= pi[i]*generator[i*n+i]; child[i] = (double)(i+1)/(double)n; }
  assert(fabsl(mean-(n > 1)) < 3e-13L);
  const double times[] = {0, 1e-12, 1e-6, 0.01, 0.5, 2, 10};
  for (unsigned mode = 0; mode < 2; ++mode) for (size_t t = 0; t < sizeof(times)/sizeof(*times); ++t) {
    if (mode == CRPX_BRANCH_FIXED) reference_exponential(q,n,times[t],reference);
    else reference_resolvent(q,n,times[t],reference);
    OK(crpx_substitution_model_transition(m,times[t],(crpx_branch_mode)mode,actual,n*n));
    matrix_invariants(actual,pi,n);
    for (size_t i = 0; i < n*n; ++i) close_value(actual[i],reference[i],2e-11L);
    OK(crpx_substitution_model_apply(m,times[t],(crpx_branch_mode)mode,child,applied,n));
    for (size_t i = 0; i < n; ++i) {
      long double value = 0;
      for (size_t j = 0; j < n; ++j) value += reference[i*n+j]*child[j];
      close_value(applied[i],value,2e-11L);
    }
  }
  free(q); free(pi); free(reference); free(actual); free(generator); free(child); free(applied);
  crpx_substitution_model_free(m);
}

static void test_binary_and_limits(void)
{
  crpx_substitution_model *m;
  OK(crpx_substitution_model_poisson(2,&m));
  double p[4];
  OK(crpx_substitution_model_transition(m,0.5,CRPX_BRANCH_FIXED,p,4));
  close_value(p[0],(1+expl(-1))/2,1e-14L);
  OK(crpx_substitution_model_transition(m,0.5,CRPX_BRANCH_SITEWISE_EXPONENTIAL,p,4));
  assert(p[0] == 0.75 && p[1] == 0.25);
  crpx_substitution_model_free(m);
  for (size_t n = 1; n <= 65; n += n < 5 ? 1 : 10) {
    OK(crpx_substitution_model_poisson(n,&m));
    double *matrix = malloc(n*n*sizeof(*matrix)), *input = malloc(n*sizeof(*input)), *output = malloc(n*sizeof(*output));
    assert(matrix && input && output);
    for (size_t i = 0; i < n; ++i) input[i] = 1;
    for (unsigned mode = 0; mode < 2; ++mode) {
      OK(crpx_substitution_model_transition(m,DBL_MAX,(crpx_branch_mode)mode,matrix,n*n));
      for (size_t i = 0; i < n*n; ++i) close_value(matrix[i],1.0L/n,1e-13L);
      OK(crpx_substitution_model_apply(m,DBL_MAX,(crpx_branch_mode)mode,input,output,n));
      for (size_t i = 0; i < n; ++i) close_value(output[i],1,1e-13L);
    }
    free(matrix); free(input); free(output); crpx_substitution_model_free(m);
  }
}

static void test_frequency_extremes(void)
{
  const double frequencies[] = {1,1e-250,1e-250,1e-250};
  for (int hky = 0; hky < 2; ++hky) {
    crpx_substitution_model *m;
    if (hky) OK(crpx_substitution_model_hky(frequencies,7,&m));
    else OK(crpx_substitution_model_f81(4,frequencies,&m));
    double q[16], p[16];
    long double pi[4], expected[16];
    reference_generator(4,frequencies,7,hky,pi,expected);
    OK(crpx_substitution_model_generator(m,q,16));
    for (size_t i = 0; i < 16; ++i) close_value(q[i],expected[i],1e-13L);
    for (unsigned mode = 0; mode < 2; ++mode) {
      OK(crpx_substitution_model_transition(m,0,(crpx_branch_mode)mode,p,16));
      for (size_t i = 0; i < 16; ++i) assert(p[i] == (i/4 == i%4));
      OK(crpx_substitution_model_transition(m,DBL_MAX,(crpx_branch_mode)mode,p,16));
      matrix_invariants(p,pi,4);
      for (size_t i = 0; i < 16; ++i) close_value(p[i],pi[i%4],1e-13L);
    }
    crpx_substitution_model_free(m);
  }
}

static void test_error_contracts(void)
{
  crpx_substitution_model *m = (crpx_substitution_model *)(uintptr_t)1;
  const double invalid[] = {1,0,1,1};
  assert(crpx_substitution_model_hky(invalid,2,&m) == CRPX_ERROR_INVALID_ARGUMENT && !m);
  assert(crpx_substitution_model_poisson(0,&m) == CRPX_ERROR_INVALID_ARGUMENT && !m);
  OK(crpx_substitution_model_poisson(2,&m));
  double matrix[] = {7,7,7,7}, child[] = {1,-1}, output[] = {7,7};
  assert(crpx_substitution_model_transition(m,-1,CRPX_BRANCH_FIXED,matrix,4) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(matrix[0] == 7 && matrix[3] == 7);
  assert(crpx_substitution_model_transition(m,1,CRPX_BRANCH_FIXED,matrix,3) == CRPX_ERROR_OUT_OF_RANGE);
  assert(matrix[0] == 7 && matrix[3] == 7);
  assert(crpx_substitution_model_apply(m,1,CRPX_BRANCH_FIXED,child,output,2) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(output[0] == 7 && output[1] == 7);
  crpx_substitution_model_free(m);
}

static void test_skewed_slow_transitions(void)
{
  /* A rare purine G can switch through an extremely short-lived pyrimidine
   * state even when its direct A/G exchangeability is tiny. In this limiting
   * generator, those two-step switches form an effective rate-1/2 process.
   * The tiny opposite-group mass must survive when the dominant mass rounds1. */
  const double frequencies[] = {1,1e-250,1e-250,1e-250};
  crpx_substitution_model *m;
  OK(crpx_substitution_model_hky(frequencies,1e-300,&m));
  double matrix[16];
  OK(crpx_substitution_model_transition(m,1,CRPX_BRANCH_FIXED,matrix,16));
  close_value(matrix[2*4],-expm1l(-0.5L),1e-13L);
  close_value(matrix[2*4+2],expl(-0.5L),1e-13L);
  OK(crpx_substitution_model_transition(m,1,CRPX_BRANCH_SITEWISE_EXPONENTIAL,matrix,16));
  close_value(matrix[2*4],1.0L/3,1e-13L);
  close_value(matrix[2*4+2],2.0L/3,1e-13L);
  OK(crpx_substitution_model_transition(m,1e-250,CRPX_BRANCH_FIXED,matrix,16));
  close_value(matrix[2*4],2e-250L*(0.25L+expm1l(-0.25L)),1e-13L);
  crpx_substitution_model_free(m);
  const double rare_purines[] = {1e-250,1,1e-250,1e-250};
  OK(crpx_substitution_model_hky(rare_purines,1e-300,&m));
  OK(crpx_substitution_model_transition(m,1e-250,CRPX_BRANCH_FIXED,matrix,16));
  close_value(matrix[2],1e-250L*(-expm1l(-0.25L)-0.25L*expl(-0.25L)),1e-13L);
  crpx_substitution_model_free(m);
}

int main(void)
{
  const double dna[] = {0.13,0.29,0.31,0.27};
  const double kappas[] = {0.01,1,2.7,1000,1e-12};
  for (size_t i = 0; i < sizeof(kappas)/sizeof(*kappas); ++i) compare_model(4,dna,kappas[i],1);
  const double binary[] = {3,7}; compare_model(2,binary,1,0);
  const double general[] = {1,2,3,7,4,8,9}; compare_model(7,general,1,0);
  const double constant[] = {9}; compare_model(1,constant,1,0);
  test_binary_and_limits(); test_frequency_extremes(); test_skewed_slow_transitions(); test_error_contracts();
  puts("analytic substitutions vs independent generators/uniformisation/resolvents: passed");
  return 0;
}
