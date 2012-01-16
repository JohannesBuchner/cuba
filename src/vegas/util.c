/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 9 Feb 05 th
*/


#include "decl.h"

static count ndim_, ncomp_;
static number neval_;
static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];
int vegasnbatch_ = 1000;
int vegasgridno_ = 0;
char vegasstate_[MAXSTATESIZE] = "";


#ifdef __GNUC__

Extern char vegasstate[MAXSTATESIZE]
  __attribute__ ((weak, alias("vegasstate_")));

Extern int vegasgridno
  __attribute__ ((weak, alias("vegasgridno_")));

Extern int vegasnbatch
  __attribute__ ((weak, alias("vegasnbatch_")));

#endif


#define SamplesAlloc(p, n) \
  MemAlloc(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))


#ifdef DEBUG
#include "debug.c"
#endif

