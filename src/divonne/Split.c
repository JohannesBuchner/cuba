/*
	Split.c
		determine optimal cuts for splitting a region
		this file is part of Divonne
		last modified 22 Jul 09 th
*/


#define BNDTOL .05
#define FRACT .5
#define BIG 1e10
#define SINGTOL 1e-4

#define LHSTOL .1
#define GAMMATOL .1

/* the next four macros must be in sync with the typedef of Bounds! */
#define Lower(d) (2*d)
#define Upper(d) (2*d + 1)
#define Dim(i) ((i) >> 1)
#define SignedDelta(i) ((i & 1) ? delta[i] : -delta[i])

typedef struct {
  count i;
  real save, delta;
  real f, df, fold;
  real lhs, row, sol;
} Cut;

typedef struct {
  real diff, err, spread;
} Errors;

typedef const Errors cErrors;

/*********************************************************************/

static inline real Div(creal a, creal b)
{
  return (b != 0 && fabs(b) < BIG*fabs(a)) ? a/b : a;
}

/*********************************************************************/

static void SomeCut(Cut *cut, Bounds *b)
{
  count dim, maxdim;
  static count nextdim = 0;
  real xmid[NDIM], ymid, maxdev;

  for( dim = 0; dim < ndim_; ++dim )
    xmid[dim] = .5*(b[dim].upper + b[dim].lower);
  ymid = Sample(xmid);

  maxdev = 0;
  maxdim = 0;
  for( dim = 0; dim < ndim_; ++dim ) {
    real ylower, yupper, dev;
    creal x = xmid[dim];
    xmid[dim] = b[dim].lower;
    ylower = Sample(xmid);
    xmid[dim] = b[dim].upper;
    yupper = Sample(xmid);
    xmid[dim] = x;

    dev = fabs(ymid - .5*(ylower + yupper));
    if( dev >= maxdev ) {
      maxdev = dev;
      maxdim = dim;
    }
  }

  if( maxdev > 0 ) nextdim = 0;
  else maxdim = nextdim++ % ndim_;

  cut->i = Upper(maxdim);
  cut->save = b[maxdim].upper;
  b[maxdim].upper = xmid[maxdim];
}

/*********************************************************************/

static inline real Volume(creal *delta)
{
  real vol = 1;
  count dim;
  for( dim = 0; dim < ndim_; ++dim )
    vol *= delta[Lower(dim)] + delta[Upper(dim)];
  return vol;
}

/*********************************************************************/

static inline real SetupEqs(Cut *cut, ccount ncut, real f)
{
  real sqsum = 0;
  Cut *c = &cut[ncut];
  while( --c >= cut ) {
    sqsum += Sq(c->lhs = f - c->f);
    f = c->f;
  }
  return sqsum;
}

/*********************************************************************/

static inline void SolveEqs(Cut *cut, count ncut,
  creal *delta, creal diff)
{
  real last = 0;
  real r = 1;
  Cut *c;

  for( c = cut; ; ++c ) {
    ccount dim = Dim(c->i);
    c->row = r -=
      Div(diff, (delta[Lower(dim)] + delta[Upper(dim)])*c->df);
    if( --ncut == 0 ) break;
    last += r*c->lhs;
  }

  last = Div(c->lhs - last, r);

  for( ; c >= cut; last += (--c)->lhs ) {
    creal delmin = -(c->delta = delta[c->i]);
    creal delmax = FRACT*(delmin + c->save);
    c->sol = Div(last, c->df);
    if( c->sol > delmax ) c->sol = .75*delmax;
    if( c->sol < delmin ) c->sol = .75*delmin;
  }
}

/*********************************************************************/

static count FindCuts(Cut *cut, Bounds *bounds, creal vol,
  real *xmajor, creal fmajor, creal fdiff)
{
  cint sign = (fdiff < 0) ? -1 : 1;

  count ncut = 0, icut;
  real delta[2*NDIM];
  real gamma, fgamma, lhssq;
  count dim, div;

  for( dim = 0; dim < ndim_; ++dim ) {
//    cBounds *b = &bounds[dim];
//    creal xsave = xmajor[dim];
    cBounds *b;
    real xsave;
    b = &bounds[dim];
    xsave = xmajor[dim];
    real dist = b->upper - xsave;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      Cut *c = &cut[ncut++];
      c->i = Upper(dim);
      c->save = dist;
      xmajor[dim] += dist *= FRACT;
      c->f = Sample(xmajor);
      xmajor[dim] = xsave;
    }
    delta[Upper(dim)] = dist;
  }

  for( dim = 0; dim < ndim_; ++dim ) {
    cBounds *b = &bounds[dim];
    creal xsave = xmajor[dim];
    real dist = xsave - b->lower;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      Cut *c = &cut[ncut++];
      c->i = Lower(dim);
      c->save = dist;
      xmajor[dim] -= dist *= FRACT;
      c->f = Sample(xmajor);
      xmajor[dim] = xsave;
    }
    delta[Lower(dim)] = dist;
  }

  if( ncut == 0 ) {
    SomeCut(cut, bounds);
    return 1;
  }

  for( ; ; ) {
    real mindiff = INFTY;
    Cut *mincut = cut;

    for( icut = 0; icut < ncut; ++icut ) {
      Cut *c = &cut[icut];
      creal diff = fabs(fmajor - c->f);
      if( diff <= mindiff ) {
        mindiff = diff;
        mincut = c;
      }
    }

    gamma = Volume(delta)/vol;
    fgamma = fmajor + (gamma - 1)*fdiff;

    if( sign*(mincut->f - fgamma) < 0 ) break;

    if( --ncut == 0 ) {
      SomeCut(cut, bounds);
      return 1;
    }

    delta[mincut->i] = mincut->save;
    memcpy(mincut, mincut + 1, (char *)&cut[ncut] - (char *)mincut);
  }

  for( icut = 0; icut < ncut; ++icut ) {
    Cut *c = &cut[icut];
    c->fold = c->f;
    c->df = (c->f - fmajor)/delta[c->i];
  }

  lhssq = SetupEqs(cut, ncut, fgamma);

repeat:
  SolveEqs(cut, ncut, delta, gamma*fdiff);

  for( div = 1; div <= 16; div *= 4 ) {
    real gammanew, lhssqnew;

    for( icut = 0; icut < ncut; ++icut ) {
      Cut *c = &cut[icut];
      real *x = &xmajor[Dim(c->i)];
      creal xsave = *x;
      delta[c->i] = c->delta + c->sol/div;
      *x += SignedDelta(c->i);
      c->f = Sample(xmajor);
      *x = xsave;
    }

    gammanew = Volume(delta)/vol;
    fgamma = fmajor + (gammanew - 1)*fdiff;
    lhssqnew = SetupEqs(cut, ncut, fgamma);

    if( lhssqnew <= lhssq ) {
      real fmax;

      if( fabs(gammanew - gamma) < GAMMATOL*gamma ) break;
      gamma = gammanew;

      fmax = fabs(fgamma);
      for( icut = 0; icut < ncut; ++icut ) {
        Cut *c = &cut[icut];
        creal dfmin = SINGTOL*c->df;
        creal sol = c->sol/div;
        real df = c->f - c->fold;
        df = (fabs(sol) < BIG*fabs(df)) ? df/sol : 1;
        c->df = (fabs(df) < fabs(dfmin)) ? dfmin : df;
        fmax = Max(fmax, fabs(c->f));
        c->fold = c->f;
      }

      if( lhssqnew < Sq((1 + fmax)*LHSTOL) ) break;
      lhssq = lhssqnew;
      goto repeat;
    }
  }

  for( icut = 0; icut < ncut; ++icut ) {
    Cut *c = &cut[icut];
    real *b = (real *)bounds + c->i;
    c->save = *b;
    *b = xmajor[Dim(c->i)] + SignedDelta(c->i);
  }

  return ncut;
}

/*********************************************************************/

static void Split(count iregion, int depth)
{
  TYPEDEFREGION;

  Cut cut[2*NDIM];
  Errors errors[NCOMP];
  count comp, ncut, nsplit, xregion, ireg, xreg;
  real tmp;

{
  Region *const region = region_ + iregion;
  selectedcomp_ = region->cutcomp;
  neval_cut_ -= neval_;
  ncut = FindCuts(cut, region->bounds, region->vol,
    (real *)region->result + region->xmajor, region->fmajor,
    region->fmajor - region->fminor);
  neval_cut_ += neval_;

  for( comp = 0; comp < ncomp_; ++comp ) {
    Errors *e = &errors[comp];
    e->diff = region->result[comp].avg;
    e->spread = e->err = 0;
  }
}

  xregion = nregions_;

  depth -= ncut;
  if( Explore(iregion, &samples_[0], depth, 1) ) {
    Cut *c;
    for( c = cut; ncut--; ++c ) {
      real *b = (real *)region_[iregion].bounds;
      ccount c0 = c->i, c1 = c0 ^ 1;
      creal tmp = b[c1];
      b[c1] = b[c0];
      b[c0] = c->save;
      if( !Explore(iregion, &samples_[0], depth++, ncut != 0) ) break;
      if( ncut ) ((real *)region_[iregion].bounds)[c1] = tmp;
    }
  }

  nsplit = nregions_ - xregion + 1;

  for( ireg = iregion, xreg = xregion; ireg < nregions_; ireg = xreg++ ) {
    cResult *result = region_[ireg].result;
    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &result[comp];
      Errors *e = &errors[comp];
      e->diff -= r->avg;
      e->err += r->err;
      e->spread += Sq(r->spread);
    }
  }

  tmp = 1./nsplit;
  for( comp = 0; comp < ncomp_; ++comp ) {
    Errors *e = &errors[comp];
    e->diff = tmp*fabs(e->diff);
    e->err = (e->err == 0) ? 1 : 1 + e->diff/e->err;
    e->spread = (e->spread == 0) ? 1 : 1 + e->diff/sqrt(e->spread);
  }

  tmp = 1 - tmp;
  for( ireg = iregion, xreg = xregion; ireg < nregions_; ireg = xreg++ ) {
    Result *result = region_[ireg].result;
    for( comp = 0; comp < ncomp_; ++comp ) {
      Result *r = &result[comp];
      cErrors *e = &errors[comp];
      creal c = tmp*e->diff;
      if( r->err > 0 ) r->err = r->err*e->err + c;
      r->spread = r->spread*e->spread + c*samples_[0].neff;
    }
  }
}

