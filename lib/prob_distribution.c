/* Deterministic probability algorithms adapted from biomcmc-lib and R 2.9.1 Rmath.
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins and the R Core Team.
 * R algorithms and their references are retained below.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "prob_distribution.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define CRPX_PD_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CRPX_PD_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define R_LN2          0.693147180559945309417232121458  /* ln(2) */
#define R_PI           3.141592653589793238462643383280  /* pi */
#define R_2PI          6.283185307179586476925286766559  /* 2*pi */
#define R_EXP_M1       0.367879441171442321595523770161  /* exp(-1) = 1/e */
#define R_SQRT_32      5.656854249492380195206754896838  /* sqrt(32) */
#define R_1_SQRT_2PI   0.398942280401432677939946059934  /* 1/sqrt(2*pi) */
#define R_LN_SQRT_2PI  0.918938533204672741780329736406  /* log(sqrt(2*pi)) = log(2*pi)/2 */
#define R_LN_SQRT_PId2 0.225791352644727432363097614947  /* log(sqrt(pi/2)) */

#define pInf INFINITY
#define mInf (-INFINITY)
#define NaN NAN
/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
static const double scalefactor = 115792089237316195423570985008687907853269984665640564039457584007913129639936.;
/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = R_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* log gamma correction factor for x >= 10 so that log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x) */
static double lgammacor (double x);
/* determines the number of terms for the double precision orthogonal series "dos" needed to insure
 * the error is no larger than "eta".  Ordinarily eta will be chosen to be one-tenth machine precision. */

/* evaluates the n-term Chebyshev series "a" at "x" */
static double chebyshev_eval (double x, const double *a, const int n);
/* This function calculates the minimum and maximum legal bounds for x in crpx_gammafn(x). These are not the only
 * bounds, but they are the only non-trivial ones to calculate. */

/* Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 * auxilary in crpx_log1pmx() and lgamma1p(); eps = relative tolerance */
static double logcf (double x, double i, double d, double eps);
/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
static double lgamma1p (double a);
/* wrapper over poisson density  (== dpois (x_P_1 - 1, lambda, log_p)) */
static double dpois_wrap (double x_plus_1, double lambda, bool log_p);
/* lowlevel poisson density (assuming E[X] = lambda, which is the scale and not rate) */
static double dpois_raw (double x, double lambda, bool log_p);
/* computes the log of the error term in Stirling's formula. */
static double stirlerr (double n);
/* Evaluates the "deviance part" bd0(x,M) =  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =  x * log(x/M) + M - x
 *  where M = E[X] = n*p (or = lambda), for   x, M > 0 */
static double bd0 (double x, double np);
/* auxiliary function for pgamma() */
static double pgamma_smallx (double x, double alph, bool log_p);
/* sum of upper terms of series (auxiliar to pgamma()) */
static double pd_upper_series (double x, double y, bool log_p);
/* sum of lower terms of series (auxiliar to pgamma()) */
static double pd_lower_series (double lambda, double y);
/* Continued fraction for calculation of (i/d) + o(i/d) ( auxiliar to pgamma())*/
static double pd_lower_cf (double i, double d);
/* Compute \f$\frac{dnorm (x, 0, 1, FALSE)}{pnorm (x, 0, 1, lower_tail = T, FALSE)} \f$ accurately */
static double dpnorm (double x, double lp);
/* Asymptotic expansion to calculate the probability that Poisson variate has value <= x. */
static double ppois_asymp (double x, double lambda, bool log_p);
/* lowlevel pgamma() calculation (without variable checking etc) */
static double pgamma_raw (double x, double alph, bool log_p);
/* chi-squared approximation to qgamma, g = log Gamma (nu/2) and tol = EPS1 */
static double qchisq_appr (double p, double nu, double g, bool log_p, double tol);
/* returns both cummulative points from the normal */
static void pnorm_both (double x, double *cum, double *ccum, int i_tail, bool log_p);
/* find quantile for poisson distribution */
static double do_poisson_search (double y, double *z, double p, double lambda, double incr);


static const double lgamma_algmcs[] = {
  +.1666389480451863247205729650822e+0,  -.1384948176067563840732986059135e-4,  +.9810825646924729426157171547487e-8,
  -.1809129475572494194263306266719e-10, +.6221098041892605227126015543416e-13, -.3399615005417721944303330599666e-15,
  +.2683181998482698748957538846666e-17, -.2868042435334643284144622399999e-19, +.3962837061046434803679306666666e-21,
  -.6831888753985766870111999999999e-23, +.1429227355942498147573333333333e-24, -.3547598158101070547199999999999e-26,
  +.1025680058010470912000000000000e-27, -.3401102254316748799999999999999e-29, +.1276642195630062933333333333333e-30
};

static const double lgamma_coeffs[] = { /* (zeta(2)-1)/2, (zeta(3)-1)/3 ... (zeta(40+1)-1)/(40+1) */
  0.3224670334241132182362075833230126e-0, 0.6735230105319809513324605383715000e-1,
  0.2058080842778454787900092413529198e-1, 0.7385551028673985266273097291406834e-2,
  0.2890510330741523285752988298486755e-2, 0.1192753911703260977113935692828109e-2,
  0.5096695247430424223356548135815582e-3, 0.2231547584535793797614188036013401e-3,
  0.9945751278180853371459589003190170e-4, 0.4492623673813314170020750240635786e-4,
  0.2050721277567069155316650397830591e-4, 0.9439488275268395903987425104415055e-5,
  0.4374866789907487804181793223952411e-5, 0.2039215753801366236781900709670839e-5,
  0.9551412130407419832857179772951265e-6, 0.4492469198764566043294290331193655e-6,
  0.2120718480555466586923135901077628e-6, 0.1004322482396809960872083050053344e-6,
  0.4769810169363980565760193417246730e-7, 0.2271109460894316491031998116062124e-7,
  0.1083865921489695409107491757968159e-7, 0.5183475041970046655121248647057669e-8,
  0.2483674543802478317185008663991718e-8, 0.1192140140586091207442548202774640e-8,
  0.5731367241678862013330194857961011e-9, 0.2759522885124233145178149692816341e-9,
  0.1330476437424448948149715720858008e-9, 0.6422964563838100022082448087644648e-10,
  0.3104424774732227276239215783404066e-10, 0.1502138408075414217093301048780668e-10,
  0.7275974480239079662504549924814047e-11, 0.3527742476575915083615072228655483e-11,
  0.1711991790559617908601084114443031e-11, 0.8315385841420284819798357793954418e-12,
  0.4042200525289440065536008957032895e-12, 0.1966475631096616490411045679010286e-12,
  0.9573630387838555763782200936508615e-13, 0.4664076026428374224576492565974577e-13,
  0.2273736960065972320633279596737272e-13, 0.1109139947083452201658320007192334e-13
};

static const double stirl_sferr_halves[] = { /* error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0. */
  0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
  0.0413406959554092940938221,  0.03316287351993628748511048, 0.02767792568499833914878929,
  0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
  0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
  0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
  0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
  0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
  0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
  0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
  0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690
};

/* log gamma correction factor for x >= 10 so that log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x) */
static double
lgammacor (double x)
{
  /* below: 94906265.62425156 = 2^26.5 and  3.745194030963158e306 = DBL_MAX / 48 =  2^1020 / 3 */
  if (x < 10) return NaN;
  if (x < 94906265.62425156) {
    double tmp = (10./x);
    return chebyshev_eval (tmp * tmp * 2 - 1, lgamma_algmcs, 5) / x;
  }
  return 1. / (x * 12.);
}

/* determines the number of terms for the double precision orthogonal series "dos" needed to insure
   the error is no larger than "eta".  Ordinarily eta will be chosen to be one-tenth machine precision. */
static double
chebyshev_eval (double x, const double *a, const int n)
{
  double b0, b1, b2, twox;
  int i;

  if ((n < 1) || (n > 1000) || (x < -1.1) || (x > 1.1)) return NaN;

  twox = x * 2;
  b2 = b1 = 0;
  b0 = 0;
  for (i = 1; i <= n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[n - i];
  }
  return (b0 - b2) * 0.5;
}

/* Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 * auxilary in crpx_log1pmx() and lgamma1p(); eps = relative tolerance */
static double
logcf (double x, double i, double d, double eps)
{
  double c1 = 2 * d, c2 = i + d, c4 = c2 + d, a1 = c2;
  double b1 = i * (c2 - i * x), b2 = d * d * x, a2 = c4 * c2 - b2;

  b2 = c4 * b1 - i * b2;

  while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
    double c3 = c2*c2*x;
    c2 += d; c4 += d;
    a1 = c4 * a2 - c3 * a1;
    b1 = c4 * b2 - c3 * b1;

    c3 = c1 * c1 * x;
    c1 += d; c4 += d;
    a2 = c4 * a1 - c3 * a2;
    b2 = c4 * b1 - c3 * b2;

    if (fabs (b2) > scalefactor) {
      a1 /= scalefactor; b1 /= scalefactor; a2 /= scalefactor; b2 /= scalefactor;
    }
    else if (fabs (b2) < 1 / scalefactor) {
      a1 *= scalefactor; b1 *= scalefactor; a2 *= scalefactor; b2 *= scalefactor;
    }
  }

  return a2 / b2;
}

static double
lgamma1p (double a)
{/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
  const double eulers_const =	 0.5772156649015328606065120900824024;
  const int N = 40; /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 */
  const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
  double lgam;
  int i;

  if (fabs (a) >= 0.5) return crpx_lgammafn (a + 1, NULL);
  /* Abramowitz & Stegun 6.1.33, also  http://functions.wolfram.com/06.11.06.0008.01 */
  lgam = c * logcf (-a / 2, N + 2, 1, 1e-14);
  for (i = N - 1; i >= 0; i--) lgam = lgamma_coeffs[i] - a * lgam;
  return (a * lgam - eulers_const) * a - crpx_log1pmx (a);
}

/* dpois_wrap (x_P_1,  lambda, log_p) == dpois (x_P_1 - 1, lambda, log_p) */
static double
dpois_wrap (double x_plus_1, double lambda, bool log_p)
{
  if (!isfinite (lambda)) return log_p ? mInf : 0.;
  if (x_plus_1 > 1.) return dpois_raw (x_plus_1 - 1, lambda, log_p);
  if (lambda > (fabs (x_plus_1 - 1) * M_cutoff))
    return log_p? (-lambda - crpx_lgammafn (x_plus_1, NULL)) : exp (-lambda - crpx_lgammafn (x_plus_1, NULL));
  else {
    double d = dpois_raw (x_plus_1, lambda, log_p);
    return log_p ? d + log (x_plus_1 / lambda) : d * (x_plus_1 / lambda);
  }
}

static double
dpois_raw (double x, double lambda, bool log_p)
{
  /* x >= 0 ; integer for dpois(), but not e.g. for pgamma()! lambda >= 0 */
  if (lambda == 0) return( (x == 0) ? log_p ? 0. : 1. : log_p ? mInf : 0. );
  if (!isfinite (lambda)) return log_p ? mInf : 0.;
  if (x < 0) return log_p ? mInf : 0.;
  if (x <= lambda * DBL_MIN) return log_p ? -lambda : exp (-lambda);
  if (lambda < x * DBL_MIN)
    return log_p ? (-lambda + x*log(lambda) - crpx_lgammafn (x+1, NULL)) :
    exp (-lambda + x*log(lambda) - crpx_lgammafn (x+1, NULL));

  return log_p ? (-0.5 * log (R_2PI * x) - stirlerr (x) - bd0 (x,lambda)) :
  exp (-stirlerr (x) - bd0 (x,lambda)) / sqrt (R_2PI * x);
}

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 * see also lgammacor() in which computes almost the same! */
static double
stirlerr (double n)
{
  static const double S[5] = {1./12. ,1./360., 1./1260., 1./1680., 1./1188};
  double nn;

  /* For n > 15, uses the series 1/12n - 1/360n^3 + ...
   * For n <=15, integers or half-integers, uses stored values.
   * For other n < 15, uses lgamma directly */
  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn) return stirl_sferr_halves[(int)nn];
    return (crpx_lgammafn (n + 1., NULL) - (n + 0.5)*log (n) + n - R_LN_SQRT_2PI);
  }
  nn = n*n;
  if (n>500.) return ((S[0]-S[1]/nn)/n);
  if (n> 80.) return ((S[0]-(S[1]-S[2]/nn)/nn)/n);
  if (n> 35.) return ((S[0]-(S[1]-(S[2]-S[3]/nn)/nn)/nn)/n);
  return ((S[0]-(S[1]-(S[2]-(S[3]-S[4]/nn)/nn)/nn)/nn)/n); /* 15 < n <= 35 : */
}

/* Evaluates the "deviance part"
 *  bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =  x * log(x/M) + M - x
 *  where M = E[X] = n*p (or = lambda), for   x, M > 0
 *  In a manner that should be stable (with small relative error) for all x and M=np. In particular for x/np close
 *  to 1, direct evaluation fails, and evaluation is based on the Taylor series of log((1+v)/(1-v)) with
 *  v = (x-np)/(x+np). */
static double
bd0 (double x, double np)
{
  double ej, s, s1, v;
  int j;

  if(!isfinite (x) || !isfinite (np) || (np == 0.)) return NaN;

  if (fabs(x-np) < 0.1*(x+np)) {
    v = (x-np)/(x+np);
    s = (x-np)*v;/* s using v -- change by MM */
    ej = 2*x*v;
    v = v*v;
    for (j=1; ; j++) { /* Taylor series */
      ej *= v;
      s1 = s+ej/((j<<1)+1);
      if (s1==s) return (s1); /* last term was effectively 0 */
      s = s1;
    }
  }
  return (x*log (x/np) + np - x); /* else:  | x - np |  is not too small */
}

static double
pgamma_smallx (double x, double alph, bool log_p)
{/* Abramowitz and Stegun 6.5.29 [right] */
  double f2, sum = 0, c = alph, n = 0, term;

  /* Relative to 6.5.29 all terms have been multiplied by alph and the first, thus being 1, is omitted. */
  do {
    n++;
    c *= -x / n;
    term = c / (alph + n);
    sum += term;
  } while (fabs (term) > DBL_EPSILON * fabs (sum));

  double f1 = log_p ? crpx_log1p (sum) : 1. + sum;
  if (alph > 1) {
    f2 = dpois_raw (alph, x, log_p);
    f2 = log_p ? f2 + x : f2 * exp (x);
  }
  else if (log_p) f2 = alph * log (x) - lgamma1p (alph);
  else f2 = pow (x, alph) / exp (lgamma1p (alph));
  return log_p ? f1 + f2 : f1 * f2;
}

static double
pd_upper_series (double x, double y, bool log_p)
{
  double term = x / y;
  double sum = term;

  /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))    =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n)) =
   *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n))) ~  x/y + o(x/y) {which happens when alph -> Inf} */
  do { y += 1.; term *= x / y; sum += term; } while (term > sum * DBL_EPSILON);
  return log_p ? log (sum) : sum;
}

static double
pd_lower_series (double lambda, double y)
{
  double term = 1, sum = 0;

  /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
   *     =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n   ~  y/lambda + o(y/lambda) */
  while (y >= 1 && term > sum * DBL_EPSILON) { term *= y / lambda; sum += term; y -= 1.; }

  if (y != floor (y)) {
    /* The series does not converge as the terms start getting bigger (besides flipping sign) for y < -lambda. */
    double f;
    f = pd_lower_cf (y, lambda + 1 - y);
    sum += term * f;
  }
  return sum;
}

static double
pd_lower_cf (double i, double d)
{/* Continued fraction for calculation of (i/d) + o(i/d) */
  static const int max_it = 200000;
  double f = 0, of, c1, c2, c3, c4, a1, b1, a2, b2;

  if (i < d * 1e-20) return (i/d); /* includes d = Inf,  or i = 0 < d */

  a1 = 0; b1 = 1;
  a2 = i; b2 = d;
  while (b2 > scalefactor) { a1/=scalefactor; b1/=scalefactor; a2/=scalefactor; b2/=scalefactor;}

  if(a2 == 0) return 0;/* just in case, e.g. d=i=0 */

  c2 = a2; c4 = b2; c1 = 0;
  while (c1 < max_it) {
    c1++;	c2--;	c3 = c1 * c2;	c4 += 2;
    a1 = c4 * a2 + c3 * a1;
    b1 = c4 * b2 + c3 * b1;

    c1++;	c2--;	c3 = c1 * c2;	c4 += 2;
    a2 = c4 * a1 + c3 * a2;
    b2 = c4 * b1 + c3 * b2;

    if (b2 > scalefactor) { a1/=scalefactor; b1/=scalefactor; a2/=scalefactor; b2/=scalefactor;}

    if (b2 != 0) {
      of = f;
      f = a2 / b2;
      /* convergence check: relative; absolute for small f : */
      if (fabs (f - of) <= DBL_EPSILON * CRPX_PD_MAX(1., fabs(f)))
        return f;
    }
  }
  return f;
}

/* Compute the following ratio with higher accuracy that would be had from doing it directly.
 *   \f$ \frac{ dnorm (x, 0, 1, FALSE) } { pnorm (x, 0, 1, FALSE) } \f$
 * Abramowitz & Stegun 26.2.12 */
static double
dpnorm (double x, double lp)
{
  /* So as not to repeat a pnorm call, we expect lp == pnorm (x, 0, 1, TRUE)
   * but use it only in the non-critical case where either x is small or p==exp(lp) is close to 1. */
  bool upper_tail = false;
  if (x < 0) { x = -x; upper_tail = true; }
  if ((x > 10) && upper_tail) {
    double term = 1. / x;
    double sum = term;
    double x2 = x * x;
    double i = 1.;

    do { term *= -i / x2; sum += term; i += 2.; } while (fabs (term) > DBL_EPSILON * sum);
    return 1. / sum;
  }
  else return crpx_dnorm (x, 0, 1, false) / exp (lp);
}

/* Asymptotic expansion to calculate the probability that Poisson variate has value <= x.
 * Various assertions about this are made (without proof) at http://members.aol.com/iandjmsmith/PoissonApprox.htm */
static double
ppois_asymp (double x, double lambda, bool log_p)
{ /* It is called by pgamma_raw with (lower_tail = false); needs change if it is possible to be called with
     lower_tail = T */
  static const double coefs_a[8] = {
    -1e99, /* placeholder used for 1-indexing */
    2/3., -4/135., 8/2835., 16/8505., -8992/12629925., -334144/492567075., 698752/1477701225. };
  static const double coefs_b[8] = {
    -1e99, /* placeholder */
    1/12., 1/288., -139/51840., -571/2488320., 163879/209018880., 5246819/75246796800., -534703531/902961561600.};
  double elfb, elfb_term;
  double res12, res1_term, res1_ig, res2_term, res2_ig;
  double dfm, pt_, s2pt, f, np;
  int i;

  dfm = lambda - x;
  /* If lambda is large, the distribution is highly concentrated about lambda.  So representation error in x or
   * lambda can lead to arbitrarily large values of pt_ and hence divergence of the coefficients of this
   * approximation. */
  pt_ = - crpx_log1pmx (dfm / x);
  s2pt = sqrt (2 * x * pt_);
  if (dfm < 0) s2pt = -s2pt;

  res12 = 0;
  res1_ig = res1_term = sqrt (x);
  res2_ig = res2_term = s2pt;
  for (i = 1; i < 8; i++) {
    res12 += res1_ig * coefs_a[i];
    res12 += res2_ig * coefs_b[i];
    res1_term *= pt_ / i ;
    res2_term *= 2 * pt_ / (2 * i + 1);
    res1_ig = res1_ig / x + res1_term;
    res2_ig = res2_ig / x + res2_term;
  }

  elfb = x;
  elfb_term = 1;
  for (i = 1; i < 8; i++) {
    elfb += elfb_term * coefs_b[i];
    elfb_term /= x;
  }
  elfb = -elfb; /* if (upper_tail) in R-mathlib but here we only call within pgamma_raw()  with upper tail) */
  f = res12 / elfb;
  np = crpx_pnorm (s2pt, 0.0, 1.0, log_p);

  if (log_p) return np + crpx_log1p (f * dpnorm (s2pt, np));
  else return np + f * crpx_dnorm (s2pt, 0., 1., log_p);
}

static double
pgamma_raw (double x, double alph, bool log_p)
{ /* Here, assume that  (x,alph) are not NA  &  alph > 0 . */
  double res, dpsum;

  if(x <= 0.)   return log_p ? mInf : 0.;
  if(x >= pInf) return log_p ?   0. : 1.;

  if (x < 1) res = pgamma_smallx (x, alph, log_p);
  else if ((x <= (alph - 1)) && (x < 0.8 * (alph + 50))) {
    /* incl. large alph compared to x */
    double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
    double d = dpois_wrap (alph, x, log_p);
    res = log_p ? sum + d : sum * d;
  }
  else if (((alph - 1) < x) && (alph < 0.8 * (x + 50))) {
    /* incl. large x compared to alph */
    double sum;
    double d = dpois_wrap (alph, x, log_p);
    if (alph < 1) {
      if (x * DBL_EPSILON > 1 - alph) sum = log_p ? 0. : 1.;
      else {
        double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
        /* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
        sum = log_p ? log (f) : f;
      }
    }
    else {
      sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
      sum = log_p ? crpx_log1p (sum) : 1 + sum;
    }
    dpsum = d + sum;
    res = log_p ? ((dpsum > -R_LN2) ? log (-expm1 (dpsum)) : crpx_log1p (-exp (dpsum))) : (1. - d * sum);
  }
  else res = ppois_asymp (alph - 1, x, log_p); /* x > 1 and x fairly near alph. */

  /* We lose a fair amount of accuracy to underflow in the cases where the final result is very close to DBL_MIN.
   * In those cases, simply redo via log space. with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
  if (!log_p && res < DBL_MIN / DBL_EPSILON) return exp (pgamma_raw (x, alph, true));
  else return res;
}

/* g = log Gamma (nu/2) and tol = EPS1 */
static double
qchisq_appr (double p, double nu, double g, bool log_p, double tol)
{
  static const double C[4] = {4.67, 6.66, 6.73, 13.32};
  double alpha, a, c, ch, p1, p2, q, t, x, rdtmacro;

  if ((isnan(p)) || (isnan(nu))) return NaN;
  if ( (log_p  && (p > 0)) || ((!log_p) && ((p < 0) || (p > 1))) ) return pInf;
  if (nu <= 0) return NaN;

  alpha = 0.5 * nu;/* = [pq]gamma() shape */
  c = alpha-1;

  if (nu < (-1.24)*(p1 = (log_p ? p : log (p)))) ch = exp ((log (alpha) + p1 + g)/alpha + R_LN2);	/*small chi-squared*/
  else if (nu > 0.32) {	/*  using Wilson and Hilferty estimate */
    x = crpx_qnorm (p, 0, 1, log_p);
    p1 = 2./(9*nu);
    ch = nu * pow (x * sqrt(p1) + 1 - p1, 3);
    /* approximation for p tending to 1: */
    rdtmacro = log_p ? (p > -R_LN2 ? log (-expm1 (p)) : crpx_log1p (-exp (p))) : (crpx_log1p(-p));
    if( ch > 2.2*nu + 6 ) ch = -2 * (rdtmacro - c * log (0.5 * ch) + g);
  }
  else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */
    ch = 0.4;
    rdtmacro = log_p ? (p > -R_LN2 ? log (-expm1 (p)) : crpx_log1p (-exp (p))) : (crpx_log1p(-p));
    a = rdtmacro + g + c*R_LN2;
    do {
      q = ch;
      p1 = 1. / (1 + ch * (C[0] + ch));
      p2 = ch * (C[2] + ch * (C[1] + ch));
      t = -0.5 + (C[0] +2 * ch) * p1 - (C[2] + ch * (C[3] + 3 * ch)) / p2;
      ch -= (1- exp (a + 0.5 * ch) * p2 * p1) / t;
    } while(fabs(q - ch) > tol * fabs(ch));
  }
  return ch;
}

void pnorm_both(double x, double *cum, double *ccum, int i_tail, bool log_p)
{
  /* i_tail in {0,1,2} means: "lower", "upper", or "both" :
     if(lower) return  *cum := P[X <= x] ; if(upper) return *ccum := P[X >  x] = 1 - P[X <= x] */
  static const double a[5] = {
    2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249, 0.065682337918207449113
  };
  static const double b[4] = {
    47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956
  };
  static const double c[9] = {
    0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226, 2494.5375852903726711,
    6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218, 1.0765576773720192317e-8
  };
  static const double d[8] = {
    22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
    18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727
  };
  static const double p[6] = {
    0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
    2.9112874951168792e-5, 0.02307344176494017303
  };
  static const double q[5] = {
    1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244, 7.29751555083966205e-5
  };
  double xden, xnum, temp, del, eps, xsq, y;
  double min = DBL_MIN;
  int i, lower, upper;

  if(isnan(x)) { *cum = *ccum = x; return; }

  eps = DBL_EPSILON * 0.5; /* Consider changing these */

  /* i_tail in {0,1,2} =^= {lower, upper, both} */
  lower = i_tail != 1; upper = i_tail != 0;

  y = fabs(x);
  if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    }
    else xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if(lower)  *cum = 0.5 + temp;
    if(upper) *ccum = 0.5 - temp;
    if(log_p) {
      if(lower)  *cum = log(*cum);
      if(upper) *ccum = log(*ccum);
    }
  }
  else if (y <= R_SQRT_32) { /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

    xsq = trunc (y * 16) / 16;
    del = (y - xsq) * (y + xsq);
    if(log_p) {
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log (temp);
      if((lower && (x > 0.)) || (upper && (x <= 0.)))
        *ccum = crpx_log1p (-exp (-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
    }
    else {
      *cum = exp (-xsq * xsq * 0.5) * exp (-del * 0.5) * temp;
      *ccum = 1.0 - *cum;
    }
    /* swap  ccum <--> cum */
    if (x > 0.) { temp = *cum; if(lower) *cum = *ccum; *ccum = temp; }
  }

  /* else	|x| > sqrt(32) = 5.657 : the next two case differentiations were really for lower=T, log=F
   * Particularly *not* for log_p ! Note that we do want symmetry(0), lower/upper -> hence use y */
  else if(log_p || (lower && -37.5193 < x  &&  x < 8.2924) || (upper && -8.2924  < x  &&  x < 37.5193) ) {
    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x);
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (R_1_SQRT_2PI - temp) / y;

    xsq = trunc (x * 16) / 16;
    del = (x - xsq) * (x + xsq);
    if(log_p) {
      *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log (temp);
      if((lower && (x > 0.)) || (upper && (x <= 0.)))
        *ccum = crpx_log1p (-exp (-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
    }
    else {
      *cum = exp (-xsq * xsq * 0.5) * exp (-del * 0.5) * temp;
      *ccum = 1.0 - *cum;
    }
    /* swap  ccum <--> cum */
    if (x > 0.) { temp = *cum; if(lower) *cum = *ccum; *ccum = temp; }
  }
  else { /* no log_p , large x such that probs are 0 or 1 */
    if(x > 0) { *cum = 1.; *ccum = 0.;	}
    else {      *cum = 0.; *ccum = 1.;	}
  }
  if(log_p) {
    if(*cum > -min)	  *cum = -0.;
    if(*ccum > -min) *ccum = -0.;
  }
  else {
    if(*cum < min)	 *cum = 0.;
    if(*ccum < min)	*ccum = 0.;
  }
  return;
}

static double
do_poisson_search (double y, double *z, double p, double lambda, double incr)
{
  double new_y;
  if(*z >= p) for(;;) { /* search to the left */
    new_y = y - incr;
    if((y == 0) || (*z = crpx_ppois (new_y, lambda, false)) < p) return y;
    y =  CRPX_PD_MAX (0, new_y);
  }
  else for(;;) { /* search to the right */
    y = y + incr;
    if ((*z = crpx_ppois (y, lambda, false)) >= p) return y;
  }
}


double
crpx_dexp_dt (double d, double lambda, double m, bool log_p)
{
  double prob_sum;
  if ( (isnan(d)) || (isnan(lambda)) || (isnan(m))|| (d < 0.) || (!isfinite(lambda)) || (!isfinite(m)) || (lambda <= 0.) || (m < 0.) || (m != floor(m))) return NaN;
  if (d > m || d < 0 || d != floor(d)) return log_p ? mInf : 0.;
  prob_sum = crpx_expm1 (-lambda * (m + 1.)) / crpx_expm1 (-lambda); /* expm1(x) = (exp(x) - 1) when x ~ 0 */

  return log_p ? ((-lambda * d) - log (prob_sum)) : (exp(-lambda * d)/prob_sum);
}

/*! \brief cdf of discrete truncated exponential (d is discrete, m is maximum value): calculates P(D <= d) */
double
crpx_pexp_dt (double d, double lambda, double m, bool log_p)
{
  if ( (isnan(d)) || (isnan(lambda)) || (isnan(m))|| (d < 0.) || (!isfinite(lambda)) || (!isfinite(m)) || (lambda <= 0.) || (m < 0.) || (m != floor(m))) return NaN;
  if (d >= m) return log_p ? 0. : 1.;
  d = floor(d);

  return log_p ?
  (log (- crpx_expm1 (-lambda * (d + 1.))) - log (- crpx_expm1 (-lambda * (m + 1.)))) :
  (crpx_expm1 (-lambda * (d + 1.)) / crpx_expm1 (-lambda * (m + 1.)));
}

/*! \brief quantile of discrete truncated exponential, that is, finds d s.t. P(D <= d) >= p */
double
crpx_qexp_dt (double p, double lambda, double m, bool log_p)
{
  double x;
  if ( (isnan(p)) || (isnan(lambda)) || (isnan(m))|| (!isfinite(lambda)) || (!isfinite(m)) || (lambda <= 0.) || (m < 0.) || (m != floor(m))) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return m; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return m;   if (p == 0.) return 0.; }
  if (log_p) p = exp (p);

  x = crpx_log1p (p * crpx_expm1 (-lambda * (m + 1.))); /* log1p(x) = log (1+x) */

  /* For tiny probabilities, subtraction can round to exactly -1. */
  return fmax(0, fmin(m, ceil(-x / lambda - 1.)));
}

double
crpx_dgamma (double x, double alpha, double beta, bool log_p)
{
  double pr;

  if ( (isnan(x)) || (isnan(alpha)) || (isnan(beta))|| (!isfinite(alpha)) || (!isfinite(beta)) || (alpha < 0.) || (beta <= 0.)) return NaN;

  if (x < 0.) return log_p ? mInf : 0.;
  if (alpha == 0.) return (x == 0.)? pInf : (log_p ? mInf : 0.);
  if (x == 0.) {
    if (alpha < 1.) return pInf;
    if (alpha > 1.) return log_p ? mInf : 0.;
    return log_p ? log(beta) : beta; /* alpha == 1. */
  }

  if (alpha < 1) {
    pr = dpois_raw (alpha, x * beta, log_p);
    return log_p ? (pr + log(alpha/x)) : (pr * alpha/x);
  }
  pr = dpois_raw (alpha-1, x * beta, log_p);
  return log_p ? (pr + log(beta)) : (pr * beta);
}

/* Compute the quantile function of the gamma distribution. This function is based on the Applied Statistics
 * Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
 * R core improvements: lower_tail, log_p; non-trivial result for p outside [0.000002, 0.999998];
 * p ~ 1 no longer gives +Inf; final Newton step(s)
 * Best, D. J. and D. E. Roberts (1975). Percentage Points of the Chi-Squared Distribution.
 * Applied Statistics 24, page 385.  */
double
crpx_qgamma (double p, double alpha, double beta, bool log_p)
{
  static const int MAXIT = 1000;
  static const double EPS1 = 1e-2, EPS2 = 5e-7, EPS_N = 1e-15; /* final precision after each step */
  static const double i420  = 1./ 420., i2520 = 1./ 2520., i5040 = 1./ 5040, pMIN = 1e-100, pMAX = (1-1e-14);
  double p_, a, b, c, g, ch, ch0, p1, p2, q, s1, s2, s3, s4, s5, s6, t, x;
  int i, max_it_Newton = 1;

  if ( (isnan(p)) || (isnan(alpha)) || (isnan(beta))|| (!isfinite(alpha)) || (!isfinite(beta)) || (alpha < 0.) || (beta <= 0.)) return NaN;

  if (log_p) { if (p > 0.) return NaN; if (!p) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (!p) return 0.; if (p==1.) return pInf; }

  if (alpha == 0) return 0.; /* all mass at 0 : */

  p_ = log_p ? exp (p) : p;
  g = crpx_lgammafn (alpha, NULL);/* log Gamma(v/2) */

  /* Phase I : Starting Approximation */
  ch = qchisq_appr (p, 2*alpha, g, log_p, EPS1); /* 2*alpha = nu = 'df' */
  if(!isfinite (ch)) { max_it_Newton = 0; goto END; }
  if(ch < EPS2) { max_it_Newton = 20; goto END; } /* Corrected according to AS 91; MM, May 25, 1999 */
  if((p_ > pMAX) || (p_ < pMIN)) { max_it_Newton = 20; goto END; } /* did return ML_POSINF or 0.;	much better: */

  /* Phase II: Iteration *	Call pgamma() [AS 239]	and calculate seven term taylor series */
  c = alpha - 1.;
  s6 = (120 + c * (346+127*c)) * i5040; /* used below, is "const" */

  ch0 = ch;/* save initial approx. */
  for (i=1; i <= MAXIT; i++) {
    q = ch;
    p1 = 0.5 * ch;
    p2 = p_ - pgamma_raw (p1, alpha, false);
    if(!isfinite (p2) || (ch <= 0)) { ch = ch0; max_it_Newton = 27; goto END; }

    t = p2 * exp (alpha*R_LN2 + g + p1 - c*log(ch));
    b = t/ch;
    a = 0.5*t - b*c;
    s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
    s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
    s3 = (210+ a*(462+a*(707+932*a))) * i2520;
    s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
    s5 = (84+2264*a + c*(1175+606*a)) * i2520;

    ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    if(fabs(q - ch) < (EPS2 * ch)) goto END;
    if(fabs(q - ch) > (0.1 * ch)) { if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;} /* diverging? also forces ch > 0 */
  }

END:
  /* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50 To: R-bugs@biostat.ku.dk
   * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
   * Improved (MM): - only if rel.Err > EPS_N (= 1e-15); optionally *iterate* Newton */
  x = 0.5*ch/beta;

  if (max_it_Newton) {
    if (!log_p) { p = log(p); log_p = true; } /* always use log scale */
    p_ = crpx_pgamma (x, alpha, beta, log_p);

    for(i = 1; i <= max_it_Newton; i++) {
      p1 = p_ - p;
      if(fabs(p1) < fabs(EPS_N * p)) break;
      if((g = crpx_dgamma (x, alpha, beta, log_p)) == mInf) break;
      /* delta x = f(x)/f'(x);
       * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
       * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x)) */
      t = x - (p1 * exp(p_ - g));/* = x - delta x" */
      p_ = crpx_pgamma (t, alpha, beta, log_p);
      if (fabs(p_ - p) > fabs(p1) || ((i > 1) && (fabs(p_ - p) == fabs(p1)))) break; /* against flip-flop */

      if(t > (1.1*x)) t = 1.1*x; /* control step length: this could have started at the initial approximation */
      else if(t < (0.9*x)) t = 0.9*x;
      x = t;
    }
  }

  return x;
}

double
crpx_pgamma (double x, double alpha, double beta, bool log_p)
{
  if ( (isnan(x)) || (isnan(alpha)) || (isnan(beta))|| (!isfinite(alpha)) || (!isfinite(beta)) || (alpha < 0.) || (beta <= 0.)) return NaN;
  x *= beta;
  if (isnan(x)) return NaN;
  if (!alpha) return (x < 0.)?(log_p? mInf: 0.) : (log_p? 0.: 1.); /* limit case; useful e.g. in pnchisq() */
  return pgamma_raw (x, alpha, log_p);
}

double
crpx_dnorm (double x, double mu, double sigma, bool log_p)
{
  if ((isnan(x)) || (isnan(mu)) || (isnan(sigma)) || ((!isfinite (x)) && (mu == x)) || (sigma < 0.)) return NaN;
  if ((!isfinite (sigma)) || ((sigma == 0.) && (x != mu))) return log_p ? mInf : 0.;
  if ((sigma == 0.) && (mu == x)) return pInf;

  x = (x - mu) / sigma;
  if(!isfinite (x)) return log_p ? mInf : 0.;

  return (log_p? -(R_LN_SQRT_2PI +  0.5 * x * x + log(sigma)) : (R_1_SQRT_2PI * exp(-0.5 * x * x)  /   sigma));
}

/*	Compute the quantile function for the normal distribution. For small to moderate probabilities, algorithm
 *	referenced below is used to obtain an initial approximation which is polished with a final Newton step.
 *	For very large arguments, an algorithm of Wichura is used.
 *	Beasley, J. D. and S. G. Springer (1977).
 *	Algorithm AS 111: The percentage points of the normal distribution, Applied Statistics, 26, 118-121.
 *  Wichura, M.J. (1988).
 *  Algorithm AS 241: The Percentage Points of the Normal Distribution. Applied Statistics, 37, 477-484. */
double
crpx_qnorm (double p, double mu, double sigma, bool log_p)
{
  double p_, q, r, val;

  if (isnan(p) || isnan(mu) || isnan(sigma) || sigma < 0 ||
      (log_p ? p > 0 : (p < 0 || p > 1))) return NaN;
  if(sigma == 0)	return mu;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return mInf; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return mInf; }

  if (log_p) p_ = exp (p);
  else p_ = p;
  q = p_ - 0.5;

  /* use AS 241: double ppnd16_(double *p, long *ifault) ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     Produces the normal deviate Z corresponding to a given lower tail area of P; Z is accurate to about 1 part
     in 10**16. */
  if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
    r = .180625 - q * q;
    val =
    q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r +
             45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
          133.14166789178437745) * r + 3.387132872796366608)
    / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r + 39307.89580009271061) * r +
           21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
        42.313330701600911252) * r + 1.);
  }
  else { /* closer than 0.075 from {0,1} boundary */
    /* r = min(p, 1-p) < 0.075 */
    if (q > 0) r = log_p ? (- crpx_expm1 (p)) : (0.5 - p + 0.5); /* 0.5 + 0.5 -p = 1 - p */
    else r = p_;

    r = sqrt (- ((log_p && (q <= 0)) ? p : log (r))); /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

    if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
      r += -1.6;
      val = (((((((r * 7.7454501427834140764e-4 + .0227238449892691845833) * r + .24178072517745061177) *
                 r + 1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) *
              r + 4.6303378461565452959) * r + 1.42343711074968357734)
      / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + .0151986665636164571966) * r +
             .14810397642748007459) * r + .68976733498510000455) * r + 1.6763848301838038494) * r +
          2.05319162663775882187) * r + 1.);
    }
    else { /* very close to  0 or 1 */
      r += -5.;
      val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + .0012426609473880784386) * r +
                 .026532189526576123093) * r + .29656057182850489123) * r + 1.7848265399172913358) * r +
              5.4637849111641143699) * r + 6.6579046435011037772)
      / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7)* r + 1.8463183175100546818e-5) * r +
             7.868691311456132591e-4) * r + .0148753612908506148525) * r + .13692988092273580531) * r +
          .59983220655588793769) * r + 1.);
    }

    if(q < 0.0) val = -val;
    /* return (q >= 0.)? r : -r ;*/
  }
  return mu + sigma * val;
}

/* The main computation evaluates near-minimax approximations derived from those in "Rational Chebyshev
 * approximations for the error function" by W. J. Cody, Math. Comp., 1969, 631-637.  This transportable program
 * uses rational functions that theoretically approximate the normal distribution function to at least 18
 * significant decimal digits.  The accuracy achieved depends on the arithmetic system, the compiler, the intrinsic
 * functions, and proper selection of the machine-dependent constants.
 * REFERENCE: Cody, W. D. (1993). ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special Function Routines
 * and Test Drivers". ACM Transactions on Mathematical Software. 19, 22-32.*/
double
crpx_pnorm (double x, double mu, double sigma, bool log_p)
{
  double p, cp;

  if ((isnan(x)) || (isnan(mu)) || (isnan(sigma)) || ((!isfinite (x)) && (mu == x)) || (sigma < 0.)) return NaN;
  if (sigma == 0.) return (x < mu) ? (log_p? mInf : 0.) : (log_p ? 0. : 1.);

  p = (x - mu) / sigma;
  if (!isfinite (p)) return (x < mu) ? (log_p? mInf : 0.) : (log_p ? 0. : 1.);

  x = p;
  pnorm_both(x, &p, &cp, 0, log_p); /* the zero stands for lower tail (0,1,2 -> lower, upper, both) */
  return p;
}

double
crpx_dlnorm (double x, double meanlog, double sdlog, bool log_p)
{
  double y;

  if ((isnan(x)) || (isnan(meanlog)) || (isnan(sdlog)) || (sdlog <= 0.)) return NaN;
  if (x <= 0.) return log_p ? mInf : 0.;

  y = (log (x) - meanlog) / sdlog;
  return (log_p ?
          -(R_LN_SQRT_2PI   + 0.5 * y * y + log (x * sdlog)) :
          R_1_SQRT_2PI * exp (-0.5 * y * y)  /  (x * sdlog));
}

double
crpx_qlnorm (double p, double meanlog, double sdlog, bool log_p)
{
  if (isnan(p) || isnan(meanlog) || isnan(sdlog) || sdlog < 0) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return 0.; }

  return exp (crpx_qnorm (p, meanlog, sdlog, log_p));
}

double
crpx_plnorm (double x, double meanlog, double sdlog, bool log_p)
{
  if ((isnan(x)) || (isnan(meanlog)) || (isnan(sdlog)) || (sdlog <= 0.)) return NaN;
  if (x <= 0.) return log_p ? mInf : 0.;
  return crpx_pnorm (log (x), meanlog, sdlog, log_p);
}

double
crpx_dpois (double x, double lambda, bool log_p)
{ /* E[X] = lambda (remember that lambda is related to "alpha" of gamma distrib, NOT "beta" */
  if ((isnan(x)) || (isnan(lambda)) || (!isfinite(lambda)) || (lambda < 0.)) return NaN;
  if ((x < 0) || (!isfinite (x)) || (fabs (x - floor (x + 0.5)) > 1e-7)) return log_p ? mInf : 0.;
  x = floor (x + 0.5);
  return dpois_raw (x, lambda, log_p);
}

/*The quantile function of the Poisson distribution. Uses the Cornish-Fisher Expansion to include a skewness
 *correction to a normal approximation.  This gives an initial value which never seems to be off by more than
 *1 or 2. A search is then conducted of values close to this initial start point [ do_poisson_search() ] . */
double
crpx_qpois (double p, double lambda, bool log_p)
{
  double mu, sigma, gamma, z, y;
  if (isnan(p) || !isfinite(lambda) || lambda < 0 || lambda > 0x1p52) return NaN;

  if (log_p) {         if (p > 0.) return NaN; if (p == 0.) return pInf; if (p == mInf) return 0.; }
  else { if ((p < 0.) || (p > 1.)) return NaN; if (p == 1.) return pInf;   if (p == 0.) return 0.; }

  if (lambda == 0) return 0.;

  mu = lambda;
  sigma = sqrt(lambda);
  gamma = 1.0/sigma; /* gamma = sigma; PR#8058 should be kurtosis which is mu^-0.5 */

  /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c - This is far from optimal [cancellation for p ~= 1, etc] */
  if (log_p) { /* need check again (cancellation!): */
    p = exp (p);
    if (p == 0.) return 0.;
    if (p == 1.) return pInf;
  }

  /* y = approx.value (Cornish-Fisher expansion) */
  z = crpx_qnorm (p, 0., 1., false);
  y = fmax(0, floor(mu + sigma * (z + gamma * (z*z - 1.) / 6.) + 0.5));
  z = crpx_ppois (y, lambda, false);
  p *= 1 - 64 * DBL_EPSILON; /* fuzz to ensure left continuity; 1 - 1e-7 may lose too much */

  /* If the mean is not too large a simple search is OK */
  if (lambda < 1e4) return do_poisson_search (y, &z, p, lambda, 1);
   {/* Otherwise be a bit cleverer in the search */
    double incr = fmax(1, floor(y * 0.001)), oldincr, new_incr;
    do {
      oldincr = incr;
      y = do_poisson_search (y, &z, p, lambda, incr);
      new_incr = floor (incr / 100.);
      incr = CRPX_PD_MAX (1, new_incr);
    } while ((oldincr > 1) && (incr > (lambda * 1e-15)));
    return y;
   }
}

double
crpx_ppois (double x, double lambda, bool log_p)
{
  if ((isnan(x)) || (isnan(lambda)) || (!isfinite(lambda)) || (lambda < 0.)) return NaN;
  if (x < 0)               return log_p ? mInf : 0.;
  if (lambda == 0.)        return log_p ? 0. : 1.;
  if (!isfinite (x)) return log_p ? 0. : 1.;
  x = floor (x + 1e-7);

  double shape = x + 1, upper;
  if (x < lambda && shape < 0.8 * (lambda + 50)) {
    /* Evaluate the upper gamma tail directly, preserving log P[X<=x] even
     * when ordinary probabilities underflow (e.g. x=0, lambda=1000). */
    upper = dpois_wrap(shape, lambda, true) + log1p(pd_lower_series(lambda, shape - 1));
  }
  else {
    double lower = crpx_pgamma(lambda, shape, 1., true);
    upper = lower > -R_LN2 ? log(-expm1(lower)) : log1p(-exp(lower));
  }
  return log_p ? upper : exp(upper);
}

double
crpx_lgammafn (double x, int *sgn)
{
  /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 */
  static const double xmax  = 2.5327372760800758e+305; /* =DBL_MAX/log(DBL_MAX)=2^1024/(1024*log(2))=2^1014/log(2)*/
  double ans, y, sinpiy;

  if (sgn != NULL) *sgn = 1;

  /* fmod (x,y) returns the floating-point remainder of dividing x by y. */
  if ((sgn != NULL) && (x < 0) && (fmod (floor (-x), 2.) == 0)) *sgn = -1;

  if (isnan(x)) return NAN;
  if (x <= 0 && x == trunc(x)) return pInf; /* Negative integer arguments are forbidden */

  y = fabs(x);
  if (y <= 10) return log (fabs (crpx_gammafn (x)));
  if (y > xmax) return pInf; /* */
  if (x > 0) {
    if(x > 1e17) return(x * (log (x) - 1.));
    else if(x > 4934720.) return(R_LN_SQRT_2PI + (x - 0.5) * log (x) - x);
    else return R_LN_SQRT_2PI + (x - 0.5) * log (x) - x + lgammacor (x); /* i.e. y = x > 10 */
  }
  /* else: x < -10; y = -x */
  sinpiy = fabs (sin (R_PI * y));
  ans = R_LN_SQRT_PId2 + (x - 0.5) * log (y) - x - log (sinpiy) - lgammacor (y);


  return ans;
}

double crpx_log1pmx (double x)
{/* Accurate calculation of log(1+x)-x, particularly for small x.  */
  static const double minLog1Value = -0.79149064;
  static const double two[4] = { 2./9., 2./7., 2./5., 2./3.};

  if (x > 1. || x < minLog1Value) return crpx_log1p(x) - x;
  else { /* expand in	 [x/(2+x)]^2 */
    double term = x / (2. + x);
    double y = term * term;
    if (fabs(x) < 1e-2) return term * ((((two[0] * y + two[1]) * y + two[2]) * y + two[3]) * y - x);
    else return term * (2 * y * logcf (y, 3, 2, 1e-14) - x);
  }
}

double crpx_gammafn(double x) { return tgamma(x); }
double crpx_log1p(double x) { return log1p(x); }
double crpx_expm1(double x) { return expm1(x); }

double crpx_logspace_add(double log_x, double log_y)
{
  if (isnan(log_x) || isnan(log_y)) return NAN;
  if (log_x == INFINITY || log_y == INFINITY) return INFINITY;
  if (log_x == -INFINITY) return log_y;
  if (log_y == -INFINITY) return log_x;
  double maximum = fmax(log_x, log_y);
  return maximum + log1p(exp(fmin(log_x, log_y) - maximum));
}

double crpx_logspace_sub(double log_x, double log_y)
{
  if (isnan(log_x) || isnan(log_y) || log_y > log_x || (log_x == INFINITY && log_y == INFINITY))
    return NAN;
  if (log_x == log_y) return -INFINITY;
  if (log_y == -INFINITY || log_x == INFINITY) return log_x;
  double difference = log_y - log_x;
  return log_x + (difference > -R_LN2 ? log(-expm1(difference)) : log1p(-exp(difference)));
}

crpx_status crpx_discrete_gamma(double alpha, double beta, double *rates, size_t count)
{
  if (!rates || !count || !isfinite(alpha) || !isfinite(beta) || alpha <= 0 || beta <= 0)
    return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  if (count == 1) {
    double mean = alpha / beta;
    if (!isfinite(mean)) return CRPX_ERROR_OVERFLOW;
    rates[0] = mean;
    return CRPX_OK;
  }
  double *means = malloc(count * sizeof(double));
  if (!means) return CRPX_ERROR_ALLOCATION;
  double previous = 0;
  for (size_t i = 0; i < count; ++i) {
    double probability = 1;
    if (i + 1 < count) {
      double boundary = crpx_qgamma((double)(i + 1) / (double)count, alpha, beta, false);
      probability = crpx_pgamma(boundary, alpha + 1, beta, false);
    }
    means[i] = (double)((long double)(probability - previous) * count * alpha / beta);
    if (!isfinite(means[i]) || means[i] < 0) { free(means); return CRPX_ERROR_OUT_OF_RANGE; }
    previous = probability;
  }
  for (size_t i = 0; i < count; ++i) rates[i] = means[i];
  free(means);
  return CRPX_OK;
}

static double uniform_positive(crpx_rng_stream *stream)
{
  double value;
  do { value = crpx_rng_stream_double(stream); } while (value == 0);
  return value;
}

double crpx_rng_norm(crpx_rng_stream *stream, double mu, double sigma)
{
  if (!stream || !isfinite(mu) || !isfinite(sigma) || sigma < 0) return NAN;
  if (sigma == 0) return mu;
  double radius = sqrt(-2 * log(uniform_positive(stream)));
  return mu + sigma * radius * cos(R_2PI * crpx_rng_stream_double(stream));
}

double crpx_rng_lnorm(crpx_rng_stream *stream, double meanlog, double sdlog)
{
  return exp(crpx_rng_norm(stream, meanlog, sdlog));
}

/* Marsaglia and Tsang, A Simple Method for Generating Gamma Variables (2000).
 * All proposal state is local; shape < 1 uses the standard power transformation. */
double crpx_rng_gamma(crpx_rng_stream *stream, double alpha, double beta)
{
  if (!stream || !isfinite(alpha) || !isfinite(beta) || alpha < 0 || beta <= 0) return NAN;
  if (alpha == 0) return 0;
  double shape = alpha < 1 ? alpha + 1 : alpha;
  double d = shape - 1.0 / 3.0, c = 1 / sqrt(9 * d);
  for (;;) {
    double normal = crpx_rng_norm(stream, 0, 1), v = 1 + c * normal;
    if (v <= 0) continue;
    v = v * v * v;
    double uniform = uniform_positive(stream), normal_squared = normal * normal;
    if (uniform < 1 - 0.0331 * normal_squared * normal_squared ||
        log(uniform) < 0.5 * normal_squared + d * (1 - v + log(v))) {
      double value = d * v / beta;
      return alpha < 1 ? value * pow(uniform_positive(stream), 1 / alpha) : value;
    }
  }
}

double crpx_rng_pois(crpx_rng_stream *stream, double lambda)
{
  /* Doubles represent every integer only up to 2^53. */
  if (!stream || !isfinite(lambda) || lambda < 0 || lambda > 0x1p52) return NAN;
  if (lambda == 0) return 0;
  return crpx_qpois(uniform_positive(stream), lambda, false);
}

struct crpx_discrete_sampler_struct {
  size_t size;
  size_t *aliases;
  double *thresholds, *probabilities;
};

void crpx_discrete_sampler_free(crpx_discrete_sampler_t sampler)
{
  if (sampler) {
    free(sampler->aliases);
    free(sampler->thresholds);
    free(sampler->probabilities);
    free(sampler);
  }
}

crpx_status crpx_discrete_sampler_new(const double *weights, size_t size, crpx_discrete_sampler_t *out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!weights || !size) return CRPX_ERROR_INVALID_ARGUMENT;
  if (size > SIZE_MAX / sizeof(double) || size > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  double maximum = 0;
  for (size_t i = 0; i < size; ++i) {
    if (!isfinite(weights[i]) || weights[i] < 0) return CRPX_ERROR_INVALID_ARGUMENT;
    maximum = fmax(maximum, weights[i]);
  }
  if (maximum == 0) return CRPX_ERROR_INVALID_ARGUMENT;
  long double total = 0;
  for (size_t i = 0; i < size; ++i) total += (long double)weights[i] / maximum;
  crpx_discrete_sampler_t sampler = calloc(1, sizeof(*sampler));
  if (!sampler) return CRPX_ERROR_ALLOCATION;
  sampler->size = size;
  sampler->aliases = malloc(size * sizeof(size_t));
  sampler->thresholds = malloc(size * sizeof(double));
  sampler->probabilities = malloc(size * sizeof(double));
  double *scaled = malloc(size * sizeof(double));
  size_t *small = malloc(size * sizeof(size_t)), *large = malloc(size * sizeof(size_t));
  if (!sampler->aliases || !sampler->thresholds || !sampler->probabilities || !scaled || !small || !large) {
    free(scaled); free(small); free(large); crpx_discrete_sampler_free(sampler);
    return CRPX_ERROR_ALLOCATION;
  }
  size_t small_count = 0, large_count = 0;
  for (size_t i = 0; i < size; ++i) {
    sampler->probabilities[i] = (double)(((long double)weights[i] / maximum) / total);
    scaled[i] = (double)((((long double)weights[i] / maximum) / total) * size);
    if (scaled[i] < 1) small[small_count++] = i;
    else large[large_count++] = i;
  }
  while (small_count && large_count) {
    size_t low = small[--small_count], high = large[--large_count];
    sampler->thresholds[low] = scaled[low];
    sampler->aliases[low] = high;
    scaled[high] = (scaled[high] - 1) + scaled[low];
    if (scaled[high] < 1) small[small_count++] = high;
    else large[large_count++] = high;
  }
  while (large_count) {
    size_t index = large[--large_count];
    sampler->thresholds[index] = 1;
    sampler->aliases[index] = index;
  }
  while (small_count) {
    size_t index = small[--small_count];
    sampler->thresholds[index] = 1;
    sampler->aliases[index] = index;
  }
  free(scaled); free(small); free(large);
  *out = sampler;
  return CRPX_OK;
}

size_t crpx_discrete_sampler_size(const struct crpx_discrete_sampler_struct *sampler)
{ return sampler ? sampler->size : 0; }

crpx_status crpx_discrete_sampler_sample(const struct crpx_discrete_sampler_struct *sampler,
                                       crpx_rng_stream *stream, size_t *out)
{
  if (!sampler || !stream || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  uint64_t column;
  crpx_status status = crpx_rng_stream_range(stream, sampler->size, &column);
  if (status != CRPX_OK) return status;
  *out = crpx_rng_stream_double(stream) < sampler->thresholds[column] ? (size_t)column : sampler->aliases[column];
  return CRPX_OK;
}

crpx_status crpx_discrete_sampler_pdf(const struct crpx_discrete_sampler_struct *sampler, size_t index, double *out)
{
  if (!sampler || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  if (index >= sampler->size) return CRPX_ERROR_OUT_OF_RANGE;
  *out = sampler->probabilities[index];
  return CRPX_OK;
}
