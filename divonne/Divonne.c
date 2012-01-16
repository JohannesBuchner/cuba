/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 13 Apr 04
*/

#include "util.c"

#define Print(s) puts(s); fflush(stdout)

static Integrand integrand_;
static PeakFinder peakfinder_;

/*********************************************************************/

static inline void DoSample(count n, ccount ldx, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, f, &phase_);
    x += ldx;
    f += ncomp_;
  }
}

/*********************************************************************/

static inline count SampleExtra(cBounds *b)
{
  count n = nextra_;
  peakfinder_(&ndim_, b, &n, xextra_);
  DoSample(n, ldxgiven_, xextra_, fextra_);
  return n;
}

/*********************************************************************/

#include "common.c"

void Divonne(ccount ndim, ccount ncomp, Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  ccount ngiven, ccount ldxgiven, real *xgiven,
  ccount nextra, PeakFinder peakfinder,
  int *pnregions, int *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( ndim < MINDIM || ndim > MAXDIM ) *pfail = -1;
  else {
    neval_ = neval_opt_ = neval_cut_ = 0;
    integrand_ = integrand;
    peakfinder_ = peakfinder;
    border_.lower = border;
    border_.upper = 1 - border_.lower;
    ngiven_ = ngiven;
    xgiven_ = NULL;
    ldxgiven_ = IMax(ldxgiven, ndim_);
    nextra_ = nextra;

    if( ngiven + nextra ) {
      ccount nxgiven = ngiven*ldxgiven;
      ccount nxextra = nextra*ldxgiven;
      ccount nfgiven = ngiven*ncomp;
      ccount nfextra = nextra*ncomp;

      Allocate(xgiven_, nxgiven + nxextra + nfgiven + nfextra);
      xextra_ = xgiven_ + nxgiven;
      fgiven_ = xextra_ + nxextra;
      fextra_ = fgiven_ + nfgiven;

      if( nxgiven ) {
        phase_ = 0;
        Copy(xgiven_, xgiven, nxgiven);
        DoSample(ngiven_, ldxgiven_, xgiven_, fgiven_);
      }
    }

    *pfail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key1, key2, key3, maxpass,
      maxchisq, mindeviation,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;

    if( xgiven_ ) free(xgiven_);
  }
}

