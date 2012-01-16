/*
	divonne-f.c
		Fortran interface for Divonne
		this file is part of Divonne
		last modified 13 Apr 04 th
*/


#include "decl.h"

#ifdef HAVE_UNDERSCORE
#define divonne divonne_
#endif

void Divonne(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  ccount ngiven, ccount ldxgiven, real *xgiven,
  ccount nextra, PeakFinder peakfinder,
  int *pnregions, int *pneval, int *pfail,
  real *integral, real *error, real *prob);


void divonne(ccount *pndim, ccount *pncomp, Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, ccount *pmineval, ccount *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  ccount *pngiven, ccount *pldxgiven, real *xgiven,
  ccount *pnextra, PeakFinder peakfinder,
  int *pnregions, int *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  Divonne(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey1, *pkey2, *pkey3, *pmaxpass,
    *pborder, *pmaxchisq, *pmindeviation,
    *pngiven, *pldxgiven, xgiven,
    *pnextra, peakfinder,
    pnregions, pneval, pfail,
    integral, error, prob);
}

