/*
	suave-f.c
		Fortran interface for Suave
		this file is part of Suave
		last modified 19 Jan 05 th
*/


#include "decl.h"

#ifdef HAVE_UNDERSCORE
#define suave suave_
#endif


void Suave(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nnew, creal flatness,
  count *pnregions, count *pneval, int *pfail,
  real *integral, real *error, real *prob);


void suave(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, ccount *pmineval, ccount *pmaxeval,
  ccount *pnnew, creal *pflatness,
  count *pnregions, count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Suave(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnnew, *pflatness,
    pnregions, pneval, pfail,
    integral, error, prob);
}

