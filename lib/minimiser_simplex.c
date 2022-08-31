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

/*! \file minimiser_simplex.c 
 *  \brief  The Simplex method of Nelder and Mead, also known as the polytope search alogorithm, from the GSL (GPL-3.0)
 * Ref: Nelder, J.A., Mead, R., Computer Journal 7 (1965) pp. 308-313. This implementation uses n+1 corner points in the simplex. 
 * This is the "simplex2" algorithm from GSL, where the size of simplex is calculated as the RMS distance of each vertex from the center 
 * rather than the mean distance. */

#include "minimiser_simplex.h"

/* Copyright (C) 2007, 2008, 2009 Brian Gough <bjg@network-theory.co.uk>, from the GSL (GPL-3.0)
 * Copyright (C) 2002 Tuomo Keskitalo <tuomo.keskitalo@iki.fi>, Ivo Alxneit  <ivo.alxneit@psi.ch>, from the GSL (GPL-3.0) */

// TODO: I assume that F() returns finite or NaN in error, but a good F() should return mInf or pInf also 
const double CRPX_pInf = 1./0., CRPX_mInf = -1./0., CRPX_NaN = 0./0.; 

typedef struct {
  double **x1;   /* simplex corner points , each row is one corner and cols are dimensions */
  double *y1,    /* function value at corner points */
         *ws1,   /* workspace 1 for algorithm */
         *ws2,   /* workspace 2 for algorithm */
         *center,/* center of all points */
         *delta, /* current step */
         *xmc;   /* x - center (workspace) */
  double S2;
  unsigned long count; /*!< number of attempts at initial state */

  /* above are from gsl_multimin_fminimizer_type (lowlevel) and below are some from gsl_multimin_fminimizer and other high level */

  double simplex_size, minimum; /*!< from gsl_multimin_fminimizer not gsl_multimin_fminimizer_type in GSL */
  double (*F)(double *, void *); /*!< function receiving n dimensions (and possibly extra parameters) and returning a double */
  void *params; /*!< extra parameters for the function */
  size_t size1, size2; /*!< size1 = n_rows (corners), size2 = n_cols (dimensions) */
  crpx_global_t cglob;
} crpx_simplex_struct, *crpx_simplex_t;

crpx_simplex_t
new_crpx_simplex_t (crpx_global_t cglob, size_t n, (*f)(double*,void*), void *extra_parameters)
{
  size_t i;
  if (!n) { crpx_logger_warning (cglob, "Cannot create a simplex of zero dimensions"); return NULL; }
  crpx_simplex_t s = (crpx_simplex_t) crpx_malloc (cglob, sizeof (crpx_simplex_struct));
  if (!s) return NULL;
  crpx_link_add_global_pointer (cglob, s->cglob); // thread-safe increase of ref_counter

  s->size1 = n+1; // number of samples (corner points)
  s->size2 = n;   // dimension of each corner point

  s->x1 = (double**) crpx_malloc (cglob, sizeof (double*) * s->size1);
  if (!s1->x1) { del_crpx_simplex_t (s); return NULL; }
  for (i = 0; i < s->size1; i++) {
    s->x1[i] = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
    if (!s->x1[i]) { del_crpx_simplex_t (s); return NULL; }
  }

  s->y1 = (double*) crpx_malloc (cglob, sizeof (double) * s->size1);
  if (!s->y1) { del_crpx_simplex_t (s); return NULL; }

  s->ws1 = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
  if (!s->ws1) { del_crpx_simplex_t (s); return NULL; }

  s->ws2 = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
  if (!s->ws2) { del_crpx_simplex_t (s); return NULL; }

  s->center = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
  if (!s->center) { del_crpx_simplex_t (s); return NULL; }

  s->delta = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
  if (!s->delta) { del_crpx_simplex_t (s); return NULL; }

  s->xmc = (double*) crpx_malloc (cglob, sizeof (double) * s->size2);
  if (!s->xmc) { del_crpx_simplex_t (s); return NULL; }

  state->count = 0;
  s->F = f;
  s->params = extra_parameters;

  return s;
}

void
del_crpx_simplex_t (crpx_simplex_t s)
{
  size_t i;
  if (!s) return;
  for (i = s->size1 - 1;; i--) crpx_free (cglob, s->x1[i]);
  crpx_free (cglob, s->x1);
  crpx_free (cglob, s->y1);
  crpx_free (cglob, s->ws1);
  crpx_free (cglob, s->ws2);
  crpx_free (cglob, s->center);
  crpx_free (cglob, s->delta);
  crpx_free (cglob, s->xmc);
  crpx_link_remove_global_pointer (s->cglob); // thread-safe decrease of ref_counter
  free (s);
}

void
local_blas_daxpy (size_t n, double a, double *x, double *y)
{
  size_t i, m = n % 4;

  for (i = 0; i < m; i++) y[i] += a * x[i];
  for (i = m; i + 3 < n; i += 4) {
    y[i] += a * x[i];
    y[i+1] += a * x[i+1];
    y[i+2] += a * x[i+2];
    y[i+3] += a * x[i+3];
  }
}
  
double
local_blas_dnrm2 (size_t n, double *x)
{
  double xabs, scale = 0.0, ssq = 1.0;
  size_t i;

  if (n == 1) return fabs (x[0]);
  for (i = 0; i < n; i++) if (x[i] != 0.0) {
    xabs = fabs (x[i]);
    if (xabs > scale) {
      ssq = 1.0 + ssq * (scale / xabs) * (scale / xabs);
      scale = xabs;
    } else {
      ssq += (xabs / scale) * (xabs / scale);
    }
  }
  return scale * sqrt (ssq);
}

static double
try_corner_move (crpx_simplex_t sim, const double coeff, size_t corner, double *xc)
{
  /* moves a simplex corner scaled by coeff (negative value represents mirroring by the middle point of the "other" corner points)
   * and gives new corner in xc and function value at xc as a return value */

  const double P = (double) sim->size1;
  /* xc = (1-coeff)*(P/(P-1)) * center(all) + ((P*coeff-1)/(P-1))*x_corner */
  double alpha = (1 - coeff) * P / (P - 1.0);
  double beta = (P * coeff - 1.0) / (P - 1.0);
  size_t i, m;

  for (size_t i = 0; i < sim->size2; i++) xc[i] = sim->center[i] * alpha; // gsl_blas_dscal (alpha,xc) with xc = center
  local_blas_daxpy (sim->size2, beta, sim->x1[corner], xc);
  return sim->F (xc, sim->params);
}

static void
update_point (crpx_simplex_t sim, size_t corner, const double *x, double val)
{
  const double xmcd = 0., P = (double) sim->size1;

  /* Compute delta = x - x_orig */
  memcpy (sim->delta, x, sizeof (double) * sim->size2); // delta is tmp array from crpx_simplex_t
  local_blas_daxpy (sim->size2, -1.0, sim->x1[corner], sim->delta);
  /* Compute xmc = x_orig - c */
  memcpy (sim->xmc, sim->x1[corner], sizeof (double) * sim->size2); // xmc is tmp array from crpx_simplex_t
  local_blas_daxpy (sim->size2, -1.0, sim->center, sim->xmc);

  /* Update size: S2' = S2 + (2/P) * (x_orig - c).delta + (P-1)*(delta/P)^2 */
  double d = local_blas_dnrm2 (sim->size2, sim->delta);
  for (size_t i = 0; i < sim->size2; i++) xmcd += sim->xmc[i] * sim->delta[i]; //  gsl_blas_ddot (state->xmc, state->delta, &xmcd);
  sim->S2 += (2.0 / P) * xmcd + ((P - 1.0) / P) * (d * d / P);

  /* Update center:  c' = c + (x - x_orig) / P */
  double alpha = 1.0 / P;
  local_blas_daxpy (-alpha, sim->x1[corner], sim->center);
  local_blas_daxpy (alpha, x, sim->center);

  memcpy (sim->x1[corner], x, sizeof (double) * sim->size2); // actual update of the corner point
  sim->y1[corner] = val; // and the function value at it
}

static bool 
contract_by_best (crpx_simplex_t sim, size_t best, double *xc)
{
/* Function contracts the simplex in respect to best valued corner. That is, all corners besides the best corner are moved.
 * (This function is rarely called in practice, since it is the lastchoice, hence not optimised - BJG)  */

  /* the xc vector is simply work space here */
  gsl_matrix *x1 = state->x1;
  gsl_vector *y1 = state->y1;

  size_t i, j;
  double status_success = true;

  for (i = 0; i < sim->size1; i++) if (i != best) {
    for (j = 0; j < x1->size2; j++) sim->x1[i][j] = 0.5 * sim->x1[i][j] + sim->x1[best][j];
    /* evaluate function in the new point */
    memcp (xc, sim->x1[i], sizeof (double) * sim->size2);
    sim->y1[i] = sim->F (xc, sim->params);
    if (sim->y1[i] == CRPX_NaN) status_success = false; /* at least one bad function value; let user handle the situation */
  }
  /* We need to update the centre and size as well */
  compute_center (sim);
  compute_size (sim);
  return status_success;
}

void
compute_center (crpx_simplex_t sim)
{  /* calculates the center of the simplex and stores in center */
  double alpha = 1.0 / (double) (sim->size1);
  size_t i;
#ifdef __STDC_IEC_559__ // almost any compiler should have this defined, it states that (int) 0 == (double) 0.0
  memset (sim->center, 0, sizeof (double) * sim->size2); // gsl_vector_set_zero (sim->center);
#else
  for (i = 0; i < sim->size2; i++) sim->center[i] = 0.0;
#endif
  for (i = 0; i < sim->size1; i++) local_blas_daxpy (sim->size2, 1.0, sim->x1[i], sim->center); // gsl_blas_daxpy (1.0, sim->x1[i], sim->center);
  for (i = 0; i < sim->size2; i++) center[i] = sim->center * alpha; // gsl_blas_dscal (alpha,center)
}

static double
compute_size (crpx_simplex_t sim) // sim->center
{ /* calculates simplex size as rms sum of length of vectors  from simplex center to corner points:     
   * sqrt( sum ( || y - y_middlepoint ||^2 ) / n )  */
  size_t i;
  double t, ss = 0.0;

  for (i = 0; i < sim->size1; i++) {
    memcpy (sim->ws1, sim->x1[i], sizeof (double) * sim->size2); // gsl_matrix_get_row (s, x1, i); with s = wc1
    local_blas_daxpy (sim->size2, -1.0, sim->center, sim->ws1);  // gsl_blas_daxpy (-1.0, sim->center, s);
    t = local_blas_dnrm2 (sim->size2, sim->ws1); // gsl_blas_dnrm2 (s);
    ss += t * t;
  }
  sim->S2 = (ss / (double) sim->size1);  /* Store squared size in the state */
  return sqrt (sim->S2);
}

size_t
crpx_simplex_initial_state (crpx_simplex_t sim, const double *x0, double *step_size)
{
  size_t i;
  double val;

  sim->count++;
  /* first point is the original x0 */
  sim->y1[0] = sim->F (x0, sim->params);
  memcpy (sim->x1[0], x0, sizeof (double) * sim->size2); // gsl_matrix_set_row (state->x1, 0, x);
  /* unlike GSL, we copy x0 even if it is NaN so that user can do a post mortem */
  if (sim->y1[0] == CRPX_NaN) return 0; // return how many successful dimensions

  /* following points are initialized to x0 + step_size */
  for (i = 0; i < sim->size2; i++) {
    memcpy (sim->ws1, x0, sizeof (double) * sim->size2); // status = gsl_vector_memcpy (xtemp, x); with xtemp=ws1
    sim->ws1[i] = x0[i] + step_size[i]; // gsl_vector_set (xtemp, i, gsl_vector_get (x, i) + gsl_vector_get (step_size, i) );
    sim->y1[i + 1] = sim->F (sim->ws1, sim->params); 
    memcpy (sim->x1[i + 1], sim->ws1, sizeof (double) * sim->size2); // gsl_matrix_set_row (state->x1, i+1, xtemp);
    if (sim->y1[i + 1] == CRPX_NaN) return i + 1; // return how many successful dimensions
  }
  compute_center (sim);
  sim->simplex_size = compute_size (sim);  /* Initialize simplex size */
  return i + 1; // sim->size2
}

//static int nmsimplex_iterate (void *vstate, gsl_multimin_function * f, gsl_vector * x, double *size, double *fval)
bool
crpx_simplex_iterate (crpx_simplex_t sim)
{ /* Simplex iteration tries to minimize function f value; w/ corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */
  /* xc and xc2 vectors store tried corner point coordinates */
  bool status;
  double val, val2;
  double dhi, ds_hi, dlo;
  size_t  hi,  s_hi,  lo;
  size_t i;

  /* get index of highest, second highest and lowest point */
  dhi = dlo = sim->y1[0];  hi = lo = 0;
  ds_hi = sim->y1[1];    s_hi = 1;

  for (i = 1; i < sim->size1; i++) {
    val = sim->y1[i];
    if (val < dlo)        {   dlo = val;   lo = i; }
    else if (val > dhi)   { ds_hi = dhi; s_hi = hi;  dhi = val;  hi = i;  }
    else if (val > ds_hi)	{ ds_hi = val; s_hi = i;  }
  }

  val = try_corner_move (sim, -1., hi, sim->ws1);  /* try reflecting the highest value point */

  if ((val != CRPX_NaN) && (val < sim->y1[lo])) { /* reflected point is lowest, try expansion */
      val2 = try_corner_move (sim, -2.0, hi, sim->ws2);
      if ((val2 != CRPX_NaN) && (val2 < sim->y1[lo])) { update_point (sim, hi, sim->ws2, val2); }
      else                                            { update_point (sim, hi, sim->ws1, val);	}

  } else if ((val == CPRX_NaN) || (val > sim->y1[s_hi])) { /* reflection does not improve things enough, or we got a non-finite function value */
    if ((val != CRPX_NaN) && (val <= sim->y1[hi])) { update_point (sim, hi, sim->ws1, val); } /* trial point is better than highest point */
    val2 = try_corner_move (sim. 0.5, hi, sim->ws2); /* try one-dimensional contraction */
    if ((val2 != CRPX_NaN) && (val2 <= sim->y1[hi])) { update_point (sim, hi, sim->ws2, val2); }
    else { /* contract the whole simplex about the best point */
      status = contract_by_best (sim, lo, sim->ws1);
      if (!status) { return false; }
    }

  } else { /* trial point is better than second highest point.  Replace highest point by it */
    update_point (sim, hi, sim->ws1, val);
  }
// STOPHERE
  /* return lowest point of simplex as x */
  lo = gsl_vector_min_index (y1); // FIXME: gsl does it by brute force 
  gsl_matrix_get_row (x, x1, lo);
  *fval = gsl_vector_get (y1, lo);

  /* Update simplex size */
  double S2 = state->S2;

  if (S2 > 0) {
    *size = sqrt (S2);
  }
  else {
    /* recompute if accumulated error has made size invalid */
    *size = compute_size (state, state->center);
  }

  return true;
}

static const gsl_multimin_fminimizer_type nmsimplex_type = 
{ "nmsimplex2",	/* name */
  sizeof (nmsimplex_state_t),
  &nmsimplex_alloc,
  &nmsimplex_set,
  &nmsimplex_iterate,
  &nmsimplex_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_nmsimplex2 = &nmsimplex_type;


static inline double
ran_unif (unsigned long *seed)
{
  unsigned long s = *seed;
  *seed = (s * 69069 + 1) & 0xffffffffUL;
  return (*seed) / 4294967296.0;
}

static int
nmsimplex_set_rand (void *vstate, gsl_multimin_function * f,
		    const gsl_vector * x,
		    double *size, const gsl_vector * step_size)
{
  size_t i, j;
  double val;

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_vector *xtemp = state->ws1;

  if (xtemp->size != x->size)
    {
      GSL_ERROR ("incompatible size of x", GSL_EINVAL);
    }

  if (xtemp->size != step_size->size)
    {
      GSL_ERROR ("incompatible size of step_size", GSL_EINVAL);
    }

  /* first point is the original x0 */

  val = GSL_MULTIMIN_FN_EVAL (f, x);

  if (!gsl_finite (val))
    {
      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
    }

  gsl_matrix_set_row (state->x1, 0, x);
  gsl_vector_set (state->y1, 0, val);

  {
    gsl_matrix_view m =
      gsl_matrix_submatrix (state->x1, 1, 0, x->size, x->size);

    /* generate a random orthornomal basis  */
    unsigned long seed = state->count ^ 0x12345678;

    ran_unif (&seed);		/* warm it up */

    gsl_matrix_set_identity (&m.matrix);

    /* start with random reflections */
    for (i = 0; i < x->size; i++)
      {
        double s = ran_unif (&seed);
        if (s > 0.5) gsl_matrix_set (&m.matrix, i, i, -1.0);
      }

    /* apply random rotations */
    for (i = 0; i < x->size; i++)
      {
	for (j = i + 1; j < x->size; j++)
	  {
	    /* rotate columns i and j by a random angle */
	    double angle = 2.0 * M_PI * ran_unif (&seed);
	    double c = cos (angle), s = sin (angle);
	    gsl_vector_view c_i = gsl_matrix_column (&m.matrix, i);
	    gsl_vector_view c_j = gsl_matrix_column (&m.matrix, j);
	    gsl_blas_drot (&c_i.vector, &c_j.vector, c, s);
	  }
      }

    /* scale the orthonormal basis by the user-supplied step_size in
       each dimension, and use as an offset from the central point x */

    for (i = 0; i < x->size; i++)
      {
	double x_i = gsl_vector_get (x, i);
	double s_i = gsl_vector_get (step_size, i);
	gsl_vector_view c_i = gsl_matrix_column (&m.matrix, i);

	for (j = 0; j < x->size; j++)
	  {
	    double x_ij = gsl_vector_get (&c_i.vector, j);
	    gsl_vector_set (&c_i.vector, j, x_i + s_i * x_ij);
	  }
      }

    /* compute the function values at each offset point */

    for (i = 0; i < x->size; i++)
      {
	gsl_vector_view r_i = gsl_matrix_row (&m.matrix, i);

	val = GSL_MULTIMIN_FN_EVAL (f, &r_i.vector);

	if (!gsl_finite (val))
	  {
	    GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
	  }

	gsl_vector_set (state->y1, i + 1, val);
      }
  }

  compute_center (state, state->center);

  /* Initialize simplex size */
  *size = compute_size (state, state->center);

  state->count++;

  return GSL_SUCCESS;
}

static const gsl_multimin_fminimizer_type nmsimplex2rand_type = 
{ "nmsimplex2rand",	/* name */
  sizeof (nmsimplex_state_t),
  &nmsimplex_alloc,
  &nmsimplex_set_rand,
  &nmsimplex_iterate,
  &nmsimplex_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_nmsimplex2rand = &nmsimplex2rand_type;

