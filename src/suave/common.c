/*
	common.c
		includes most of the modules
		this file is part of Suave
		last modified 28 Apr 13 th
*/


#define RegionAlloc(t, p, n, nnew) { \
  const size_t size = sizeof(Region) + \
    (n)*SAMPLESIZE + \
    (nnew)*t->ndim*sizeof(bin_t); \
  MemAlloc(p, size); \
  p->size = size; \
}

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > NDIM ) return true;
  return t->ndim < SOBOL_MINDIM ||
    (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"
#include "Sample.c"
#include "Fluct.c"

