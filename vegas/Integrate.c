/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 26 Nov 04
*/


static int Integrate(creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nstart, ccount nincrease,
  real *integral, real *error, real *prob)
{
  real *sample;
  count dim, comp;
  int fail = 1;
  struct {
    count niter, nsamples, neval;
    Cumulants cumul[NCOMP];
    Grid grid[NDIM];
  } state;
  int statemsg = VERBOSE;
  struct stat st;

  if( VERBOSE > 1 ) {
    char s[256];
    sprintf(s, "Vegas input parameters:\n"
      "  ndim %d\n  ncomp %d\n"
      "  epsrel %g\n  epsabs %g\n"
      "  flags %d\n  mineval %d\n  maxeval %d\n"
      "  nstart %d\n  nincrease %d",
      ndim_, ncomp_,
      epsrel, epsabs,
      flags, mineval, maxeval,
      nstart, nincrease);
    Print(s);
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  IniRandom(2*maxeval, ndim_);

  if( *vegasstate_ && stat(vegasstate_, &st) == 0 &&
      st.st_size == sizeof(state) && (st.st_mode & 0400) ) {
    cint h = open(vegasstate_, O_RDONLY);
    read(h, &state, sizeof(state));
    close(h);
    SkipRandom(neval_ = state.neval);

    if( VERBOSE ) {
      char s[256];
      sprintf(s, "\nRestoring state from %s.", vegasstate_);
      Print(s);
    }
  }
  else {
    state.niter = 0;
    state.nsamples = nstart;
    Zap(state.cumul);
    GetGrid(state.grid);
  }

  /* main iteration loop */

  for( ; ; ) {
    count n;
    creal jacobian = 1./state.nsamples;
    real *w, *x, *f, *lastf;
    bin_t *bins;

    SamplesAlloc(sample, state.nsamples);
    w = sample;
    x = w + state.nsamples;
    f = x + state.nsamples*ndim_;
    lastf = f + state.nsamples*ncomp_;
    bins = (bin_t *)lastf;

    for( n = state.nsamples; n; --n ) {
      real weight = jacobian;

      GetRandom(x);

      for( dim = 0; dim < ndim_; ++dim ) {
        creal pos = *x*NBINS;
        ccount bin = pos;
        creal prev = (bin == 0) ? 0 : state.grid[dim][bin - 1];
        creal diff = state.grid[dim][bin] - prev; 
        *x++ = prev + (pos - bin)*diff;
        *bins++ = bin;
        weight *= diff*NBINS;
      }

      *w++ = weight;
    }

    DoSample(state.nsamples, w, f);

    w = sample;

    while( f < lastf ) {
      creal weight = *w++;

      for( comp = 0; comp < ncomp_; ++comp ) {
        Cumulants *c = &state.cumul[comp];
        creal wfun = weight*(*f++);
        c->sum += wfun;
        c->sqsum += Sq(wfun);
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &state.cumul[comp];
      real avg, sigsq;
      real w = Weight(c->sum, c->sqsum, state.nsamples);

      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > Max(epsabs, fabs(c->avg)*epsrel));

      if( state.niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration %d:  %d integrand evaluations so far",
        state.niter + 1, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cCumulants *c = &state.cumul[comp];
        p += sprintf(p, "\n[%d] %g +- %g  \tchisq %g (%d df)",
          comp + 1, c->avg, c->err, c->chisq, state.niter);
      }

      Print(s);
    }

    if( fail == 0 && neval_ >= mineval ) {
      if( *vegasstate_ ) unlink(vegasstate_);
      break;
    }

    if( neval_ >= maxeval && *vegasstate_ == 0 ) break;

    Reweight(state.grid, sample, x, f, state.cumul);

    ++state.niter;
    state.nsamples += nincrease;

    if( *vegasstate_ ) {
      cint h = creat(vegasstate_, 0666);
      if( h != -1 ) {
        state.neval = neval_;
        write(h, &state, sizeof(state));
        close(h);

        if( statemsg ) {
          char s[256];
          sprintf(s, "\nSaving state to %s.", vegasstate_);
          Print(s);
          statemsg = false;
        }
      }
    }

    if( neval_ >= maxeval ) break;

    free(sample);
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cCumulants *c = &state.cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, state.niter);
  }

abort:
  free(sample);
  PutGrid(state.grid);

  return fail;
}

