/*
	decl.h
		Type declarations
		this file is part of Vegas
		last modified 18 Nov 04 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>


/* to compile with a non-gnu cc, set the following to fixed values */
#ifndef NDIM
#define NDIM ndim_
#endif
#ifndef NCOMP
#define NCOMP ncomp_
#endif


#define VERBOSE (flags & 3)
#define LAST (flags & 4)

#define NOTZERO 0x1p-104

#define NBINS 128

#define MAXGRIDS 10

#define STATESIZE 128


typedef enum { false, true } bool;

typedef const bool cbool;

typedef unsigned char bin_t;

typedef const bin_t cbin_t;

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

typedef real Grid[NBINS];

typedef struct {
  real sum, sqsum;
  real weightsum, avgsum;
  real chisum, chisqsum, guess;
  real avg, err, chisq;
} Cumulants;

typedef const Cumulants cCumulants;


typedef const void (*Integrand)(ccount *, creal *, ccount *, real *);

