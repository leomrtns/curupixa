/* Private shared log-space transition kernels. SPDX-License-Identifier: MIT */
/*! \file likelihood_kernel_internal.h
    \brief Internal transition and projection calculations shared by likelihood and ancestral reconstruction.
 */
#ifndef CURUPIXA_LIKELIHOOD_KERNEL_INTERNAL_H
#define CURUPIXA_LIKELIHOOD_KERNEL_INTERNAL_H
#include "substitution_model_internal.h"
#include "likelihood.h"
typedef struct {
  double loss, factor, matrix[16];
} crpx_log_kernel;
double crpx_private_likelihood_log_add (double a, double b);
void crpx_private_likelihood_prepare_kernel (const crpx_substitution_model *model, double parameter,
                                             crpx_branch_mode mode, crpx_log_kernel *kernel);
void crpx_private_likelihood_project (const crpx_substitution_model *model, const crpx_log_kernel *kernel,
                                      const double *child, const double *log_pi, double *out);
/* Analytic models use reversibility; general rates use a direct transpose. scratch and out each contain states entries and must not overlap input or each other.
   Computes P^T times a parent message. */
void crpx_private_likelihood_project_transpose (const crpx_substitution_model *model, const crpx_log_kernel *kernel,
                                                const double *parent, const double *log_pi, double *scratch,
                                                double *out);
struct crpx_likelihood_workspace;
const crpx_log_kernel *crpx_private_likelihood_kernels (const struct crpx_likelihood_workspace *workspace);
double crpx_private_likelihood_log_transition (const crpx_substitution_model *model, const crpx_log_kernel *kernel,
                                               size_t parent, size_t child);
void crpx_private_likelihood_project_max (const crpx_substitution_model *model, const crpx_log_kernel *kernel,
                                          const double *child, double *out);
crpx_status crpx_private_likelihood_evaluate_joint (const crpx_tree *, const crpx_character_data *,
                                                    const crpx_substitution_model *, const crpx_likelihood_options *,
                                                    crpx_likelihood_workspace *, double *);
#endif
