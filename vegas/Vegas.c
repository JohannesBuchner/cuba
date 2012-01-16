/*
	Vegas.c
		Vegas Monte-Carlo integration
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

void Vegas(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nstart, ccount nincrease, 
  count *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( ndim < MINDIM || ndim > MAXDIM ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate(epsrel, epsabs,
      flags, mineval, maxeval, nstart, nincrease,
      integral, error, prob);

    *pneval = neval_;
  }
}

