/*
	CSample.c
		the serial sampling routine
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 11 Apr 14 th
*/


workerini cubaini;


static inline number SampleRaw(cThis *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  number nvec = t->nvec;
  for( ; n > 0; n -= nvec ) {
    nvec = IMin(n, nvec);
    if( t->integrand(&t->ndim, x, &t->ncomp, f, t->userdata, &nvec
          VES_ONLY(, w, &iter)
          DIV_ONLY(, &t->phase)) == ABORT ) return -1;
    VES_ONLY(w += nvec;)
    x += nvec*t->ndim;
    f += nvec*t->ncomp;
  }
  return 0;
}

/*********************************************************************/

static inline void DoSampleSerial(This *t, cnumber n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  if( t->initfun ) {
    t->initfun(cubaini.initarg);
    t->initfun = NULL;
    t->exitfun = cubaini.exitfun;
  }
  t->neval += n;
  if( SampleRaw(t, n, x, f VES_ONLY(, w, iter)) ) 
    longjmp(t->abort, -99);
}

/*********************************************************************/

#ifdef HAVE_FORK

static void DoSample(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter));
DIV_ONLY(static int Explore(This *t, cint iregion);)

#else

#define DoSample DoSampleSerial
#define Explore ExploreSerial
#define ForkCores(t)
#define WaitCores(t)

#endif

#ifdef DIVONNE
static inline count SampleExtra(This *t, cBounds *b)
{
  number n = t->nextra;
  t->peakfinder(&t->ndim, b, &n, t->xextra);
  DoSample(t, n, t->xextra, t->fextra);
  return n;
}
#endif

#include "common.c"

#ifdef HAVE_FORK
#include "Fork.c"
#endif

#include "Integrate.c"

