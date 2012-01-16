/*
	util.c
		Utility functions
		this file is part of Cuhre
		last modified 17 Jan 05 th
*/


#include "decl.h"

static count ndim_, ncomp_, nregions_;
static number neval_;


#define Allocate(p, n) \
  if( (p = malloc((n)*sizeof(*p))) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }


#ifdef DEBUG
#include "debug.c"
#endif

