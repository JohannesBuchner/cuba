/*
	util.c
		Utility functions
		this file is part of Suave
		last modified 17 Jan 05 th
*/


#include "decl.h"

static count ndim_, ncomp_, nregions_;
static number neval_;


#define Allocate(p, n) \
  if( (p = malloc(n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define RegionAlloc(p, n, nnew) \
  Allocate(p, sizeof(Region) + \
              (n)*(ndim_ + ncomp_ + 1)*sizeof(real) + \
              (nnew)*ndim_*sizeof(bin_t))


#ifdef DEBUG
#include "debug.c"
#endif

