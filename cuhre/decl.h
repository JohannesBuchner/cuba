/*
	decl.h
		Type declarations
		this file is part of Cuhre
		last modified 13 Apr 04 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


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


typedef enum { false, true } bool;

typedef const bool cbool;

typedef unsigned char bin_t;

typedef const bin_t cbin_t;

typedef /*unsigned*/ int count;

typedef const count ccount;

typedef const int cint;

typedef /*long*/ double real;
	/* Switching to long double is not as trivial as it
	   might seem here.  sqrt, erf, exp, pow need to be
	   replaced by their long double versions (sqrtl, ...),
	   printf formats need to be updated similarly, and
	   ferrying long doubles to Mathematica is of course
	   quite another matter, too. */

typedef const real creal;


typedef struct {
  real avg, err;
  count bisectdim;
} Result;

typedef const Result cResult;


typedef struct {
  real avg, err, lastavg, lasterr;
  real weightsum, avgsum;
  real guess, chisum, chisqsum, chisq;
} Totals;

typedef const Totals cTotals;


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real *x, *f;
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;


#define TYPEDEFREGION \
  typedef struct region { \
    struct region *next; \
    count div; \
    Result result[NCOMP]; \
    Bounds bounds[NDIM]; \
  } Region


typedef const void (*Integrand)(ccount *, creal *, ccount *, real *);

