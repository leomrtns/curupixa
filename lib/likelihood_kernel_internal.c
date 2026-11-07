/* Shared by pruning and full-data reconstruction. SPDX-License-Identifier: GPL-3.0-or-later */
#include "likelihood_kernel_internal.h"
#include <math.h>
double crpx_private_likelihood_log_add(double a, double b)
{
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  if (a < b) { double swap=a; a=b; b=swap; }
  return a+log1p(exp(b-a));
}
void crpx_private_likelihood_prepare_kernel(const crpx_substitution_model *m, double t,
                                            crpx_branch_mode mode, crpx_log_kernel *k)
{
  if (m->hky) crpx_private_hky_log_matrix(m,t,mode,k->matrix);
  else if (!t || !m->base_decay) { k->factor=0; k->loss=-INFINITY; }
  else {
    double logx=log(m->base_decay)+log(t);
    if (mode == CRPX_BRANCH_SITEWISE_EXPONENTIAL) {
      double denom=crpx_private_likelihood_log_add(0,logx);
      k->factor=-denom; k->loss=logx-denom;
    } else {
      double x=m->base_decay*t;
      k->factor=-x; k->loss=x?log(-expm1(-x)):logx;
    }
  }
}
void crpx_private_likelihood_project(const crpx_substitution_model *m, const crpx_log_kernel *k,
                                     const double *child, const double *log_pi, double *out)
{
  if (m->hky) {
    for (size_t i=0;i<4;++i) {
      double value=-INFINITY;
      for (size_t j=0;j<4;++j) value=crpx_private_likelihood_log_add(value,k->matrix[i*4+j]+child[j]);
      out[i]=value;
    }
  } else {
    double value=-INFINITY;
    for (size_t j=0;j<m->states;++j) value=crpx_private_likelihood_log_add(value,log_pi[j]+child[j]);
    for (size_t i=0;i<m->states;++i)
      out[i]=crpx_private_likelihood_log_add(k->loss+value,k->factor+child[i]);
  }
}
void crpx_private_likelihood_project_transpose(const crpx_substitution_model *m, const crpx_log_kernel *k,
                                               const double *parent, const double *log_pi,
                                               double *scratch, double *out)
{
  for (size_t i=0;i<m->states;++i) scratch[i]=parent[i]-log_pi[i];
  crpx_private_likelihood_project(m,k,scratch,log_pi,out);
  for (size_t i=0;i<m->states;++i) out[i]+=log_pi[i];
}
