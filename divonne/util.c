/*
	util.c
		Utility functions
		this file is part of Divonne
		last modified 13 Apr 04 th
*/


#include "decl.h"

static count ndim_, ncomp_, selectedcomp_, nregions_;
static int neval_, neval_opt_, neval_cut_, sign_, phase_;

static Bounds border_;

static Samples samples_[3];
static Rule rule7_, rule9_, rule11_, rule13_;
static real *xgiven_, *fgiven_, *xextra_, *fextra_;
static count ngiven_, nextra_, ldxgiven_;

static Totals *totals_;


#define Allocate(p, n) \
  if( (p = malloc((n)*sizeof(*p))) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define VecCopy(d, s) memcpy(d, s, ndim_*sizeof(*(d)))

#define ResCopy(d, s) memcpy(d, s, ncomp_*sizeof(*(d)))

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define ResClear(d) memset(d, 0, ncomp_*sizeof(*(d)))

#define Zap(d) memset(d, 0, sizeof(d))


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


#define BITS (sizeof(int)*8 - 1)

static inline int NegQ(cint a)
{
/* return (a < 0) ? -1 : 0; */
  return a >> BITS;
}

static inline int IDim(cint a)
{
/* return (a < 0) ? 0 : a; */
  return a & NegQ(-a);
}

static inline int IMin(cint a, cint b)
{
/* return (a < b) ? a : b; */
  return a - IDim(a - b);
}

static inline int IMax(cint a, cint b)
{
/* return (a > b) ? a : b; */
  return b + IDim(a - b);
}

static inline int TrueQ(cint a)
{
/* return (a == 0) ? 0 : -1; */
  return NegQ(a | -a);
}

static inline int Min1(cint a)
{
/* return a + (a == 0); */
  return a + 1 + TrueQ(a);
}

static inline int Abs1(cint a)
{
/* return fabs(a) + (a == 0); */
  return (a ^ NegQ(a)) - NegQ(a - 1);
}


#ifdef DEBUG
#include "debug.c"
#endif

