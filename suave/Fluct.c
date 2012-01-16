/*
	Fluct.c
		compute the fluctuation in the left and right half
		this file is part of Suave
		last modified 29 Apr 04
*/


typedef struct {
  long double fluct;
  count n;
} Var;

/*********************************************************************/

static void Fluct(Var *var, real flatness,
  cBounds *b, creal *w, count n, ccount comp, creal avg, creal err)
{
  creal *x = w + n, *f = x + n*ndim_ + comp;
  creal max = ldexp(1, (LDBL_MAX_EXP - 2)/flatness);
  creal norm = 1/(err*Max(fabs(avg), err));
  count nvar = 2*ndim_;

  Clear(var, nvar);

  while( n-- ) {
    count dim;
    const long double t =
      powl(Min(1 + fabs(*w++)*Sq(*f - avg)*norm, max), flatness);

    f += ncomp_;

    for( dim = 0; dim < ndim_; ++dim ) {
      Var *v = &var[2*dim + (*x++ >= b[dim].mid)];
      const long double f = v->fluct + t;
      v->fluct = (f > LDBL_MAX/2) ? LDBL_MAX/2 : f;
      ++v->n;
    }
  }

  flatness = 2/3./flatness;
  while( nvar-- ) {
    var->fluct = powl(var->fluct, flatness);
    ++var;
  }
}

