/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 20 Jan 05 th
*/


#include "decl.h"

static count ndim_, ncomp_;
static number neval_;
static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];
int vegasgridno_ = 0;
char vegasstate_[MAXSTATESIZE] = "";


#ifdef __GNUC__

extern char vegasstate[MAXSTATESIZE]
  __attribute__ ((weak, alias("vegasstate_")));

extern int vegasgridno
  __attribute__ ((weak, alias("vegasgridno_")));

#endif


#define Allocate(p, n) \
  if( (p = malloc(n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define SamplesAlloc(p, n) \
  Allocate(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))


#ifdef DEBUG
#include "debug.c"
#endif

