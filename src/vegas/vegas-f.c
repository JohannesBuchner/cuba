/*
	vegas-f.c
		Fortran interface for Vegas
		this file is part of Vegas
		last modified 4 Jan 06 th
*/


#include "decl.h"

Extern char vegasstate_[MAXSTATESIZE];

#ifdef HAVE_UNDERSCORE
#define vegas vegas_
#endif


Extern void Vegas(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  number *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void vegas(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnstart, cnumber *pnincrease, 
  number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  /* make sure the filename is null-terminated */
  if( *vegasstate_ ) {
    char *p;
    vegasstate_[sizeof(vegasstate_) - 1] = 0;
    if( (p = strchr(vegasstate_, ' ')) ) *p = 0;
  }

  Vegas(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnstart, *pnincrease,
    pneval, pfail,
    integral, error, prob);
}

