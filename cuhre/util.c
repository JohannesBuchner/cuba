/*
	util.c
		Utility functions
		this file is part of Cuhre
		last modified 13 Apr 04 th
*/


#include "decl.h"

static count ndim_, ncomp_, nregions_, neval_;


#define Allocate(p, n) \
  if( (p = malloc((n)*sizeof(*p))) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define Zap(d) memset(d, 0, sizeof(d))


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

static inline int IMax(cint a, cint b)
{
/* return (a > b) ? a : b; */
  return b + IDim(a - b);
}


static inline real Max(creal a, creal b)
{
  return (a > b) ? a : b;
}

static inline real Sq(creal x)
{
  return x*x;
}


#ifdef DEBUG
#include "debug.c"
#endif

