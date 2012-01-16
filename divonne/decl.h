/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 3 Jun 04 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>


/* to compile with a non-gnu cc, set the following to fixed values */
#ifndef NDIM
#define NDIM ndim_
#endif
#ifndef NCOMP
#define NCOMP ncomp_
#endif


#define VERBOSE (flags & 3)
#define LAST (flags & 4)
#define REGIONS (flags & 256)

#define INFTY HUGE_VAL

#define NOTZERO 0x1p-104

#define EXTRAPOLATE_EPS (.25*border_.lower)
/*#define EXTRAPOLATE_EPS 0x1p-26*/


typedef enum { false, true } bool;

typedef const bool cbool;

typedef /*unsigned*/ int count;

typedef const count ccount;

typedef const int cint;

typedef const long clong;

typedef /*long*/ double real;
	/* Switching to long double is not as trivial as it
	   might seem here.  sqrt, erf, exp, pow need to be
	   replaced by their long double versions (sqrtl, ...),
	   printf formats need to be updated similarly, and
	   ferrying long doubles to Mathematica is of course
	   quite another matter, too. */

typedef const real creal;


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real avg, spread, maxspread, secondspread;
  void *maxregion;
} Totals;


typedef struct {
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;


typedef struct samples {
  real weight;
  real *x, *f, *avg, *err;
  void (*sampler)(const struct samples *, cBounds *, creal);
  cRule *rule;
  count coeff;
  int n, neff;
} Samples;

typedef const Samples cSamples;


#define TYPEDEFREGION \
  typedef struct { \
    real avg, err, spread, chisq; \
    real fmin, fmax; \
    real xmin[NDIM], xmax[NDIM]; \
  } Result; \
  typedef const Result cResult; \
  typedef struct region { \
    struct region *next; \
    count cutcomp, depth; \
    real *xmajor, fmajor, fminor, vol; \
    Bounds bounds[NDIM]; \
    Result result[NCOMP]; \
  } Region


typedef const void (*Integrand)(ccount *, creal *, ccount *, real *, cint *);

typedef const void (*PeakFinder)(ccount *, cBounds *, count *, real *);

