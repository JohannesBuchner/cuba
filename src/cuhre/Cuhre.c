/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
		last modified 1 May 13 th
*/


#define CUHRE
#define ROUTINE "Cuhre"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key, cchar *statefile,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  t.ndim = ndim;
  t.ncomp = ncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = flags;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.key = key;
  t.statefile = statefile;
 
  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}

/*********************************************************************/

Extern void EXPORT(cuhre)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey, cchar *statefile,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob, cint statefilelen)
{
  This t;
  t.ndim = *pndim;
  t.ncomp = *pncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.epsrel = *pepsrel;
  t.epsabs = *pepsabs;
  t.flags = *pflags;
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.key = *pkey;
  t.statefile = CString(statefile, statefilelen);

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}
