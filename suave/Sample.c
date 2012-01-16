/*
	Sample.c
		the sampling step of Suave
		this file is part of Suave
		last modified 9 Apr 04
*/


typedef struct {
  real sum, sqsum;
  real weight, weightsum, avg, avgsum;
  real guess, chisum, chisqsum;
} Cumulants;

/*********************************************************************/

static void Sample(ccount nnew, void *voidregion,
  real *lastw, real *lastx, real *lastf, cint flags)
{
  TYPEDEFREGION;

  Region *const region = voidregion;
  count comp, dim, n, df;
  Cumulants cumul[NCOMP];
  char **ss, *s;
  ccount chars = 128*(region->div + 1);

  creal jacobian = 1/ldexp(nnew, region->div);
  real *w = region->w, *f = lastx;
  bin_t *bins = (bin_t *)(lastf + nnew*ncomp_);

  for( n = nnew; n; --n ) {
    real weight = jacobian;

    GetRandom(f);

    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      creal pos = *f*NBINS;
      ccount bin = pos;
      creal prev = (bin == 0) ? 0 : b->grid[bin - 1];
      creal diff = b->grid[bin] - prev;
      *f++ = b->lower + (prev + (pos - bin)*diff)*(b->upper - b->lower);
      *bins++ = bin;
      weight *= diff*NBINS;
    }

    *lastw++ = weight;
  }

  DoSample(nnew, lastx, lastf);

  *(lastw - 1) = -*(lastw - 1);
  region->n = lastw - w;

  if( VERBOSE > 2 ) {
    char *p0;
    Allocate(ss, ndim_*64 + ncomp_*(sizeof(char *) + chars));
    s = (char *)(ss + ncomp_);
    p0 = s + ndim_*64;
    for( comp = 0; comp < ncomp_; ++comp ) {
      ss[comp] = p0;
      p0 += chars;
    }
  }

  Zap(cumul);
  df = n = 0;

  while( w < lastw ) {
    cbool final = (*w < 0);
    creal weight = fabs(*w++);
    ++n;

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &cumul[comp];

      creal wfun = weight*(*f++);
      c->sum += wfun;
      c->sqsum += Sq(wfun);

      if( final ) {
        if( n > 1 ) {
          real w = Weight(c->sum, c->sqsum, n);
          c->weightsum += c->weight = w;
          c->avgsum += c->avg = w*c->sum;

          if( VERBOSE > 2 ) {
            creal sig = sqrt(1/w);
            ss[comp] += (df == 0) ?
              sprintf(ss[comp], "\n[%d] %g +- %g (%d)", comp + 1,
                c->sum, sig, n) :
              sprintf(ss[comp], "\n    %g +- %g (%d)",
                c->sum, sig, n);
          }

          if( df == 0 ) c->guess = c->sum;
          else {
            c->chisum += w *= c->sum - c->guess;
            c->chisqsum += w*c->sum;
          }
        }

        c->sum = c->sqsum = 0;
      }
    }

    if( final ) ++df, n = 0;
  }

  region->df = --df;

  for( comp = 0; comp < ncomp_; ++comp ) {
    Result *r = &region->result[comp];
    Cumulants *c = &cumul[comp];
    creal sigsq = 1/c->weightsum;
    creal avg = sigsq*c->avgsum;

    if( LAST ) {
      r->sigsq = 1/c->weight;
      r->avg = r->sigsq*c->avg;
    }
    else {
      r->sigsq = sigsq;
      r->avg = avg;
    }
    r->err = sqrt(r->sigsq);

    r->chisq = (sigsq < .9*NOTZERO) ? 0 : c->chisqsum - avg*c->chisum;
      /* This catches the special case where the integrand is constant
         over the entire region.  Unless that constant is zero, only the
         first set of samples will have zero variance, and hence weight
         (n - 1) 1e30 (see above).  All other sets have been sampled
         from a non-constant weight function and therefore inevitably
         show some variance.  This is an artificial effect, brought about
         by the fact that the constancy of the integrand in the region is
         seen only in this subdivision, and can degrade the chi-square
         score quite a bit.  If the constancy was determined from more
         than two samples (hence .9*NOTZERO), the chi-squares from the
         other sets are removed here. */
  }

  if( VERBOSE > 2 ) {
    char *p = s;
    char *p0 = p + ndim_*64;

    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      p += sprintf(p, 
        (dim == 0) ? "\nRegion (%f) - (%f)" :
                     "\n       (%f) - (%f)",
        b->lower, b->upper);
    }

    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &region->result[comp];
      p += sprintf(p, "%s  \tchisq %g (%d df)", p0, r->chisq, df);
      p0 += chars;
    }

    Print(s);
    free(ss);
  }
}

