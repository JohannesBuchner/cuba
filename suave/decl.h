/*
	decl.h
		Type declarations
		this file is part of Suave
		last modified 16 Jul 04 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>


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

#define NBINS 64
#define MINSAMPLES 10


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

typedef real Grid[NBINS];

typedef const Grid cGrid;

typedef struct {
  real avg, err, sigsq, chisq;
} Result;

typedef const Result cResult;


typedef struct {
  real lower, upper, mid;
  Grid grid;
} Bounds;

typedef const Bounds cBounds;


#define TYPEDEFREGION \
  typedef struct region { \
    struct region *next; \
    count div, n, df; \
    Result result[NCOMP]; \
    Bounds bounds[NDIM]; \
    real fluct[NCOMP][NDIM][2]; \
    real w[0]; \
  } Region


typedef const void (*Integrand)(ccount *, creal *, ccount *, real *);

