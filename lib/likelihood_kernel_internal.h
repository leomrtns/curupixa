/* Private shared log-space transition kernels. SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_LIKELIHOOD_KERNEL_INTERNAL_H
#define CURUPIXA_LIKELIHOOD_KERNEL_INTERNAL_H
#include "substitution_model_internal.h"
typedef struct { double loss, factor, matrix[16]; } crpx_log_kernel;
double crpx_private_likelihood_log_add(double a, double b);
void crpx_private_likelihood_prepare_kernel(const crpx_substitution_model *model, double parameter,
                                            crpx_branch_mode mode, crpx_log_kernel *kernel);
void crpx_private_likelihood_project(const crpx_substitution_model *model, const crpx_log_kernel *kernel,
                                     const double *child, const double *log_pi, double *out);
/* All models are reversible. scratch and out each contain states entries and
 * must not overlap input or each other. Computes P^T times a parent message. */
void crpx_private_likelihood_project_transpose(const crpx_substitution_model *model,
                                               const crpx_log_kernel *kernel, const double *parent,
                                               const double *log_pi, double *scratch, double *out);
#endif
