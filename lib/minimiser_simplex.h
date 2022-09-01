/* This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file minimiser_simplex.h 
 *  \brief  The Simplex method of Nelder and Mead, also known as the polytope search alogorithm, from the GSL (GPL-3.0)
 * Ref: Nelder, J.A., Mead, R., Computer Journal 7 (1965) pp. 308-313. This implementation uses n+1 corner points in the simplex. 
 * This is the "simplex2" algorithm from GSL, where the size of simplex is calculated as the RMS distance of each vertex from the center 
 * rather than the mean distance. */

#ifndef _curupixa_minimiser_simplex_h_
#define _curupixa_minimiser_simplex_h_
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "global/global_variable.h"

/* from the Gnu Scientific Library (GPL-3.0): 
 *   Copyright (C) 2007, 2008, 2009 Brian Gough <bjg@network-theory.co.uk>
 *   Copyright (C) 2002 Tuomo Keskitalo <tuomo.keskitalo@iki.fi>, Ivo Alxneit  <ivo.alxneit@psi.ch> */

typedef struct {
  double **x1;   /* simplex corner points , each row is one corner and cols are dimensions */
  double *y1,    /* function value at corner points */
         *ws1,   /* workspace 1 for algorithm */
         *ws2,   /* workspace 2 for algorithm */
         *center,/* center of all points */
         *delta, /* current step */
         *xmc;   /* x - center (workspace) */
  double S2; /*!< simplex_size squared */
  unsigned long count; /*!< number of attempts at initial state */
  /* above are from gsl_multimin_fminimizer_type (lowlevel) and below are some from gsl_multimin_fminimizer and other high level */
  double *min_x, simplex_size, min_y; /*!< from gsl_multimin_fminimizer not gsl_multimin_fminimizer_type in GSL */
  double (*F)(double *, void *); /*!< function receiving n dimensions (and possibly extra parameters) and returning a double */
  void *params; /*!< extra parameters for the function */
  size_t size1, size2; /*!< size1 = n_rows (corners), size2 = n_cols (dimensions) */
  crpx_global_t cglob;
} crpx_simplex_struct, *crpx_simplex_t;

/*! \brief  Create a new simplex structure for function minimisation.
 *  \param  n_rows number of rows in the simplex
 *  \param  n_cols number of columns in the simplex
 *  \param  F      function to minimize, receiving vector of n dimensions (and possibly extra parameters), returning a double (CRPX_NaN if error)
 *  \param  extra_parameters extra parameters for the function (e.g. number of dimensions)
 *  \return a new simplex structure */
crpx_simplex_t new_crpx_simplex_t (crpx_global_t cglob, size_t n, (*f)(double*,void*), void *extra_parameters);
/*! \brief  Free a simplex structure.
 *  \param  sim simplex structure to free */
void del_crpx_simplex_t (crpx_simplex_t sim);
/*! \brief  Initialise a simplex structure creating n+1 points (for n dimensions) where each point differ from initial state in one dimension.
 *  \param  sim       simplex structure to initialise
 *  \param  x         initial guess for the minimiser
 *  \param  step_size vector with step size (epslon) for each dimension 
 *  \return the number of successfully calculated vertices (points); should equal n+1 */
size_t crpx_simplex_initial_state (crpx_simplex_t sim, const double *x0, double *step_size);
/*! \brief update parameters of currently initialised simplex structure
  *  \param  sim simplex structure to update
  *  \param  extra_parameters extra parameters for the function (e.g. number of dimensions) 
  *  \return the number of successfully calculated vertices (points); should equal n+1 */
size_t crpx_simplex_update_params (crpx_simplex_t sim, void *extra_parameters);
/*! \brief  Run the simplex algorithm for one iteration (updates one vertex of the polytope per iteration).
 *  \param  sim simplex structure to run 
 *  \return true if update was successful; the `sim->simplex_size` variable should decrease, meaning that points get closer together */
bool crpx_simplex_iterate (crpx_simplex_t sim);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
