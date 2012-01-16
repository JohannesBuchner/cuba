/*
	cuhre-f.c
		Fortran interface for Cuhre
		this file is part of Cuhre
		last modified 4 Jan 06 th
*/


#include "decl.h"

#ifdef HAVE_UNDERSCORE
#define cuhre cuhre_
#endif


Extern void Cuhre(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob);


Extern void cuhre(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Cuhre(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
    integral, error, prob);
}

