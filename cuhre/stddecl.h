/*
	stddecl.h
		Type declarations common to all Cuba routines
		last modified 17 Jan 05 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
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
#define PSEUDORNG (flags & 8)
#define REGIONS (flags & 256)

#define INFTY HUGE_VAL

#define NOTZERO 0x1p-104


#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define VecCopy(d, s) Copy(d, s, ndim_)

#define ResCopy(d, s) Copy(d, s, ncomp_)

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define VecClear(d) Clear(d, ndim_)

#define ResClear(d) Clear(d, ncomp_)

#define Zap(d) memset(d, 0, sizeof(d))


typedef enum { false, true } bool;

typedef const bool cbool;

typedef const int cint;

typedef const long clong;

#define COUNT "%d"
typedef /*unsigned*/ int count;
typedef const count ccount;

#define NUMBER "%d"
#define NUMBER7 "%7d"
typedef /*long long*/ int number;
typedef const number cnumber;

#define REAL "%g"
#define REALF "%f"
typedef /*long*/ double real;
	/* Switching to long double is not as trivial as it
	   might seem here.  sqrt, erf, exp, pow need to be
	   replaced by their long double versions (sqrtl, ...),
	   printf formats need to be updated similarly, and
	   ferrying long doubles to Mathematica is of course
	   quite another matter, too. */

typedef const real creal;


static inline real Sq(creal x)
{
  return x*x;
}

static inline real Min(creal a, creal b)
{
  return (a < b) ? a : b;
}

static inline real Max(creal a, creal b)
{
  return (a > b) ? a : b;
}

static inline real Weight(creal sum, creal sqsum, cnumber n)
{
  creal w = sqrt(sqsum*n);
  return (n - 1)/Max((w + sum)*(w - sum), NOTZERO);
}


/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))

/* (a < 0) ? 0 : a */
#define IDim(a) ((a) & NegQ(-(a)))

/* (a < b) ? a : b */
#define IMin(a, b) ((a) - IDim((a) - (b)))

/* (a > b) ? a : b */
#define IMax(a, b) ((b) + IDim((a) - (b)))

/* (a == 0) ? 0 : -1 */
#define TrueQ(a) NegQ((a) | (-a))

/* a + (a == 0) */
#define Min1(a) ((a) + 1 + TrueQ(a))

/* abs(a) + (a == 0) */
#define Abs1(a) (((a) ^ NegQ(a)) - NegQ((a) - 1))

