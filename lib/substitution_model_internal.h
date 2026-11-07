/* Private closed-form kernel cache. SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_SUBSTITUTION_MODEL_INTERNAL_H
#define CURUPIXA_SUBSTITUTION_MODEL_INTERNAL_H
#include "substitution_model.h"
#include <stdbool.h>
struct crpx_substitution_model {
  size_t states;
  bool hky;
  double *pi, base_decay, transition_decay, group_decay[2], group_mass[2], kappa;
};
void crpx_private_hky_log_matrix(const crpx_substitution_model *, double, crpx_branch_mode, double [16]);
#endif
