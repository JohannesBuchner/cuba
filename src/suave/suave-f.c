/*
	suave-f.c
		Fortran interface for Suave
		this file is part of Suave
		last modified 4 Jan 06 th
*/


#include "decl.h"

#ifdef HAVE_UNDERSCORE
#define suave suave_
#endif


Extern void Suave(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void suave(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnnew, creal *pflatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Suave(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnnew, *pflatness,
    pnregions, pneval, pfail,
    integral, error, prob);
}

