/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 18 Nov 04 th
*/


#include "decl.h"

static count ndim_, ncomp_, neval_;
static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];
char vegasstate_[STATESIZE];
int vegasgridno_;


#define Allocate(p, n) \
  if( (p = malloc(n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define SamplesAlloc(p, n) \
  Allocate(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))

#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define Zap(d) memset(d, 0, sizeof(d))


static inline real Max(creal a, creal b)
{
  return (a > b) ? a : b;
}

static inline real Sq(creal x)
{
  return x*x;
}

static inline real Weight(creal sum, creal sqsum, ccount n)
{
  creal w = sqrt(sqsum*n);
  return (n - 1)/Max((w + sum)*(w - sum), NOTZERO);
}

#ifdef DEBUG
#include "debug.c"
#endif

