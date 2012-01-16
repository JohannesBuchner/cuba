/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 25 May 09 th
*/


#include "stddecl.h"

#define EXTRAPOLATE_EPS (.25*border_.lower)
/*#define EXTRAPOLATE_EPS 0x1p-26*/


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real avg, spreadsq;
  real spread, secondspread;
  real nneed, maxerrsq, mindevsq;
  int iregion;
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
  number n, neff;
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
    count cutcomp, depth, xmajor; \
    real fmajor, fminor, vol; \
    Bounds bounds[NDIM]; \
    Result result[NCOMP]; \
  } Region

#define CHUNKSIZE 4096


typedef void (*Integrand)(ccount *, creal *, ccount *, real *, cint *);

typedef void (*PeakFinder)(ccount *, cBounds *, number *, real *);

