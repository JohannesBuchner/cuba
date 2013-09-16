/*
	Explore.c
		sample region, determine min and max, split if necessary
		this file is part of Divonne
		last modified 2 Aug 13 th
*/


typedef struct {
  real fmin, fmax;
  creal *xmin, *xmax;
} Extrema;

/*********************************************************************/

static int ExploreSerial(This *t, ccount iregion)
{
  csize_t regionsize = RegionSize;
  Region *region = RegionPtr(iregion);
  cBounds *bounds = region->bounds;
  Result *result = RegionResult(region);

  Vector(Extrema, extrema, NCOMP);
  Vector(real, xtmp, NDIM);
  Result *r, *r0;
  creal *x;
  real *f;
  real halfvol, maxerr;
  count n, dim, comp, maxcomp;
  cSamples *samples = &t->samples[region->isamples];

  for( comp = 0; comp < t->ncomp; ++comp ) {
    Extrema *e = &extrema[comp];
    e->fmin = INFTY;
    e->fmax = -INFTY;
    e->xmin = e->xmax = NULL;
  }

  if( region->isamples == 0 ) {		/* others already sampled */
    real vol = 1;
    for( dim = 0; dim < t->ndim; ++dim ) {
      cBounds *b = &bounds[dim];
      vol *= b->upper - b->lower;
    }
    region->vol = vol;

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Result *r = &result[comp];
      r->fmin = INFTY;
      r->fmax = -INFTY;
    }

    x = t->xgiven;
    f = t->fgiven;
    n = t->ngiven;
    if( t->nextra ) n += SampleExtra(t, bounds);

    for( ; n; --n ) {
      for( dim = 0; dim < t->ndim; ++dim ) {
        cBounds *b = &bounds[dim];
        if( x[dim] < b->lower || x[dim] > b->upper ) goto skip;
      }
      for( comp = 0; comp < t->ncomp; ++comp ) {
        Extrema *e = &extrema[comp];
        creal y = f[comp];
        if( y < e->fmin ) e->fmin = y, e->xmin = x;
        if( y > e->fmax ) e->fmax = y, e->xmax = x;
      }
skip:
      x += t->ldxgiven;
      f += t->ncomp;
    }

    samples->sampler(t, iregion);
  }

  x = samples->x;
  f = samples->f;
  for( n = samples->n; n; --n ) {
    for( comp = 0; comp < t->ncomp; ++comp ) {
      Extrema *e = &extrema[comp];
      creal y = *f++;
      if( y < e->fmin ) e->fmin = y, e->xmin = x;
      if( y > e->fmax ) e->fmax = y, e->xmax = x;
    }
    x += t->ndim;
  }
  t->neval_opt -= t->neval;

  halfvol = .5*region->vol;
  maxerr = -INFTY;
  maxcomp = -1;

  for( comp = 0; comp < t->ncomp; ++comp ) {
    Extrema *e = &extrema[comp];
    Result *r = &result[comp];
    real ftmp, err;

    if( e->xmin ) {	/* not all NaNs */
      t->selectedcomp = comp;
      XCopy(xtmp, e->xmin);
      ftmp = FindMinimum(t, bounds, xtmp, e->fmin);
      if( ftmp < r->fmin ) {
        r->fmin = ftmp;
        XCopy(&r->xminmax[0], xtmp);
      }

      t->selectedcomp = Tag(comp);
      XCopy(xtmp, e->xmax);
      ftmp = -FindMinimum(t, bounds, xtmp, -e->fmax);
      if( ftmp > r->fmax ) {
        r->fmax = ftmp;
        XCopy(&r->xminmax[t->ndim], xtmp);
      }
    }

    r->spread = halfvol*(r->fmax - r->fmin);
    err = r->spread/Max(fabs(r->avg), NOTZERO);
    if( err > maxerr ) {
      maxerr = err;
      maxcomp = comp;
    }
  }

  t->neval_opt += t->neval;

  if( maxcomp == -1 ) { /* all NaNs */
    region->depth = 0;
    return -1;
  }

  region->cutcomp = maxcomp;
  r0 = RegionResult(region);
  r = r0 + maxcomp;
  if( halfvol*(r->fmin + r->fmax) > r->avg ) {
    region->fminor = r->fmin;
    region->fmajor = r->fmax;
    region->xmajor = &r->xminmax[t->ndim] - (real *)r0;
  }
  else {
    region->fminor = r->fmax;
    region->fmajor = r->fmin;
    region->xmajor = &r->xminmax[0] - (real *)r0;
  }

  if( region->isamples == 0 ) {
    if( (region->depth < INIDEPTH && r->spread < samples->neff*r->err) ||
        r->spread < t->totals[maxcomp].secondspread ) region->depth = 0;
    if( region->depth == 0 )
      for( comp = 0; comp < t->ncomp; ++comp )
        t->totals[comp].secondspread =
          Max(t->totals[comp].secondspread, result[comp].spread);
  }

  if( region->depth ) Split(t, iregion);

  return iregion;
}

