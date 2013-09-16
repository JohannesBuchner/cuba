/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Suave
		checkpointing by B. Chokoufe
		last modified 2 May 13 th
*/


static int Integrate(This *t, real *integral, real *error, real *prob)
{
  TYPEDEFREGION;

  count dim, comp;
  Region *anchor = NULL, *region = NULL;
  int fail;
  struct {
    signature_t signature;
    count nregions, df;
    number neval;
    Result totals[NCOMP];
  } state;
  StateDecl;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Suave input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nnew " NUMBER "\n  flatness " REAL "\n"
      "  statefile \"%s\"\n",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->nnew, t->flatness,
      t->statefile);
    Print(s);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  ShmAlloc(t, ShmRm(t));
  ForkCores(t);

  if( (fail = setjmp(t->abort)) ) goto abort;

  t->epsabs = Max(t->epsabs, NOTZERO);
  IniRandom(t);

  StateSetup(t);

  if( StateReadTest(t) ) {
    StateReadOpen(t, fd) {
      count iregion;
      size_t *regsize, totsize;
      Region **last = &anchor;
      if( read(fd, &state, sizeof state) != sizeof state ||
        state.signature != StateSignature(t, 2) ) break;
      t->nregions = state.nregions;
      totsize = t->nregions*sizeof *regsize;
      regsize = alloca(totsize);
      if( read(fd, regsize, totsize) != totsize ) break;
      for( iregion = 0; iregion < t->nregions; ++iregion )
        totsize += regsize[iregion];
      if( st.st_size != sizeof state + totsize ) break;
      for( iregion = 0; iregion < t->nregions; ++iregion ) {
        Region *reg;
        MemAlloc(reg, regsize[iregion]);
        read(fd, reg, regsize[iregion]);
        *last = reg;
        last = &reg->next;
      }
      *last = NULL;
    } StateReadClose(t, fd);
    t->neval = state.neval;
    t->rng.skiprandom(t, t->neval);
  }

  if( ini ) {
    RegionAlloc(t, anchor, t->nnew, t->nnew);
    anchor->next = NULL;
    anchor->div = 0;
    t->nregions = 1;

    for( dim = 0; dim < t->ndim; ++dim ) {
      Bounds *b = &anchor->bounds[dim];
      b->lower = 0;
      b->upper = 1;
      b->mid = .5;

      if( dim == 0 ) {
        count bin;
        /* define the initial distribution of bins */
        for( bin = 0; bin < NBINS; ++bin )
          b->grid[bin] = (bin + 1)/(real)NBINS;
      }
      else Copy(b->grid, anchor->bounds[0].grid, NBINS);
    }

    t->neval = 0;
    Sample(t, t->nnew, anchor, anchor->w,
      anchor->w + t->nnew,
      anchor->w + t->nnew + t->ndim*t->nnew);
    state.df = anchor->df;
    FCopy(state.totals, anchor->result);
  }

  /* main iteration loop */
  for( ; ; ) {
    Var var[NDIM][2], *vLR;
    real maxratio, maxerr, minfluct, bias, mid;
    Region *regionL, *regionR, *reg, **parent, **par;
    Bounds *bounds, *boundsL, *boundsR;
    count maxcomp, bisectdim;
    number n, nL, nR, nnewL, nnewR;
    real *w, *wL, *wR, *x, *xL, *xR, *f, *fL, *fR, *wlast, *flast;

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        t->nregions, t->neval);

      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *tot = &state.totals[comp];
        p += sprintf(p, "\n[" COUNT "] " 
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, tot->avg, tot->err, tot->chisq, state.df);
      }

      Print(s);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( comp = 0; comp < t->ncomp; ++comp ) {
      creal ratio = state.totals[comp].err/MaxErr(state.totals[comp].avg);
      if( ratio > maxratio ) {
        maxratio = ratio;
        maxcomp = comp;
      }
    }

    if( maxratio <= 1 && t->neval >= t->mineval ) break;

    if( t->neval >= t->maxeval ) {
      fail = 1;
      break;
    }

    maxerr = -INFTY;
    parent = &anchor;
    region = anchor;
    for( par = &anchor; (reg = *par); par = &reg->next ) {
      creal err = reg->result[maxcomp].err;
      if( err > maxerr ) {
        maxerr = err;
        parent = par;
        region = reg;
      }
    }

    /* find the bisectdim which minimizes the fluctuations */
    Fluct(t, var[0],
      region->bounds, region->w, region->n, maxcomp,
      region->result[maxcomp].avg, Max(maxerr, t->epsabs));

    bias = (t->epsrel < 1e-50) ? 2 :
      Max(pow(2., -(real)region->div/t->ndim)/t->epsrel, 2.);
    minfluct = INFTY;
    bisectdim = 0;
    for( dim = 0; dim < t->ndim; ++dim ) {
      cBounds *b = &region->bounds[dim];
      creal fluct = (var[dim][0].fluct + var[dim][1].fluct)*
        (bias - b->upper + b->lower);
      if( fluct < minfluct ) {
        minfluct = fluct;
        bisectdim = dim;
      }
    }

    /* apply stratified sampling to distribute points in bisected region */
    vLR = var[bisectdim];
    minfluct = vLR[0].fluct + vLR[1].fluct;
    nnewL = IMax(
      (minfluct == 0) ? t->nnew/2 : (count)(vLR[0].fluct/minfluct*t->nnew),
      MINSAMPLES );
    nL = vLR[0].n + nnewL;
    nnewR = IMax(t->nnew - nnewL, MINSAMPLES);
    nR = vLR[1].n + nnewR;

    RegionAlloc(t, regionL, nL, nnewL);
    RegionAlloc(t, regionR, nR, nnewR);

    *parent = regionL;
    regionL->next = regionR;
    regionR->next = region->next;
    regionL->div = regionR->div = region->div + 1;

    bounds = &region->bounds[bisectdim];
    mid = bounds->mid;
    n = region->n;
    w = wlast = region->w;  x = w + n;     f = flast = x + n*t->ndim;
    wL = regionL->w;        xL = wL + nL;  fL = xL + nL*t->ndim;
    wR = regionR->w;        xR = wR + nR;  fR = xR + nR*t->ndim;

    while( n-- ) {
      cbool final = (*w < 0);
      if( x[bisectdim] < mid ) {
        if( final && wR > regionR->w ) wR[-1] = -fabs(wR[-1]);
        *wL++ = *w++;
        XCopy(xL, x);
        xL += t->ndim;
        FCopy(fL, f);
        fL += t->ncomp;
      }
      else {
        if( final && wL > regionL->w ) wL[-1] = -fabs(wL[-1]);
        *wR++ = *w++;
        XCopy(xR, x);
        xR += t->ndim;
        FCopy(fR, f);
        fR += t->ncomp;
      }
      x += t->ndim;
      f += t->ncomp;
      if( n && final ) wlast = w, flast = f;
    }

    Reweight(t, region->bounds, wlast, flast, f, state.totals);
    XCopy(regionL->bounds, region->bounds);
    XCopy(regionR->bounds, region->bounds);

    boundsL = &regionL->bounds[bisectdim];
    boundsR = &regionR->bounds[bisectdim];
    boundsL->mid = .5*(boundsL->lower + (boundsL->upper = mid));
    boundsR->mid = .5*((boundsR->lower = mid) + boundsR->upper);

    StretchGrid(bounds->grid, boundsL->grid, boundsR->grid);

    Sample(t, nnewL, regionL, wL, xL, fL);
    Sample(t, nnewR, regionR, wR, xR, fR);

    state.df += regionL->df + regionR->df - region->df;

    for( comp = 0; comp < t->ncomp; ++comp ) {
      cResult *r = &region->result[comp];
      Result *rL = &regionL->result[comp];
      Result *rR = &regionR->result[comp];
      Result *tot = &state.totals[comp];
      real diff, sigsq;

      tot->avg += diff = rL->avg + rR->avg - r->avg;

      diff = Sq(.25*diff);
      sigsq = rL->sigsq + rR->sigsq;
      if( sigsq > 0 ) {
        creal c = Sq(1 + sqrt(diff/sigsq));
        rL->sigsq *= c;
        rR->sigsq *= c;
      }
      rL->err = sqrt(rL->sigsq += diff);
      rR->err = sqrt(rR->sigsq += diff);

      tot->sigsq += rL->sigsq + rR->sigsq - r->sigsq;
      tot->err = sqrt(tot->sigsq);

      tot->chisq += rL->chisq + rR->chisq - r->chisq;
    }

    free(region);
    region = NULL;
    ++t->nregions;

    if( StateWriteTest(t) ) {
      StateWriteOpen(t, fd) {
        Region *reg;
        size_t regsize[t->nregions], *s = regsize;
        state.signature = StateSignature(t, 2);
        state.neval = t->neval;
        state.nregions = t->nregions;
        for( reg = anchor; reg; reg = reg->next )
          *s++ = reg->size;
        write(fd, &state, sizeof state);
        write(fd, regsize, sizeof regsize);
        for( reg = anchor; reg; reg = reg->next )
          write(fd, reg, reg->size);
      } StateWriteClose(t, fd);
    }
  }

  for( comp = 0; comp < t->ncomp; ++comp ) {
    cResult *tot = &state.totals[comp];
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, state.df);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", t->nregions);
    for( region = anchor; region; region = region->next ) {
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < t->ndim; ++dim ) {
        cBounds *b = &region->bounds[dim];
        lower[dim] = b->lower;
        upper[dim] = b->upper;
      }

      MLPutFunction(stdlink, "Cuba`Suave`region", 4);
      MLPutRealList(stdlink, lower, t->ndim);
      MLPutRealList(stdlink, upper, t->ndim);

      MLPutFunction(stdlink, "List", t->ncomp);
      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->err, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->df);
    }
  }
#endif

abort:
  free(region);
  while( (region = anchor) ) {
    anchor = anchor->next;
    free(region);
  }
  WaitCores(t);
  ShmFree(t);

  StateRemove(t);

  return fail;
}

