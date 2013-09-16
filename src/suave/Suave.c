/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 2 May 13 th
*/


#define SUAVE
#define ROUTINE "Suave"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Suave)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness,
  cchar *statefile,
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
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.nnew = nnew;
  t.flatness = flatness;
  t.statefile = statefile;

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}

/*********************************************************************/

Extern void EXPORT(suave)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnnew, creal *pflatness,
  cchar *statefile,
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
  t.seed = *pseed;
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.nnew = *pnnew;
  t.flatness = *pflatness;
  t.statefile = CString(statefile, statefilelen);

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}

