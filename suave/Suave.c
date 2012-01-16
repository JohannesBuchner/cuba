/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 13 Apr 04
*/


#include "util.c"

#define Print(s) puts(s); fflush(stdout)

static Integrand integrand_;

/*********************************************************************/

static inline void DoSample(count n, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, f);
    x += ndim_;
    f += ncomp_;
  }
}

/*********************************************************************/

#include "common.c"

void Suave(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nnew, creal flatness,
  count *pnregions, count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( ndim < MINDIM || ndim > MAXDIM ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, nnew, flatness,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

