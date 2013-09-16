/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Cuhre
		checkpointing by B. Chokoufe
		last modified 2 May 13 th
*/


#define POOLSIZE 1024

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  TYPEDEFREGION;
  typedef struct pool {
    struct pool *next;
    Region region[POOLSIZE];
  } Pool;

  count dim, comp, ipool, npool;
  Pool *cur = NULL, *pool;
  ccount poolsize = sizeof(Pool);
  Region *region;
  int fail;
  struct {
    signature_t signature;
    count nregions, ncur;
    number neval;
    Totals totals[NCOMP];
  } state;
  StateDecl;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Cuhre input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key " COUNT "\n"
      "  statefile \"%s\"",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->mineval, t->maxeval,
      t->key,
      t->statefile);
    Print(s);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  t->epsabs = Max(t->epsabs, NOTZERO);

  RuleAlloc(t);
  t->mineval = IMax(t->mineval, t->rule.n + 1);
  FrameAlloc(t, ShmRm(t));
  ForkCores(t);

  if( (fail = setjmp(t->abort)) ) goto abort;

  StateSetup(t);

  if( StateReadTest(t) ) {
    StateReadOpen(t, fd) {
      Pool *prev = NULL;
      int size;
      if( read(fd, &state, sizeof state) != sizeof state ||
        state.signature != StateSignature(t, 4) ) break;
      t->neval = state.neval;
      t->nregions = state.nregions;
      do {
        Alloc(cur, 1);
        cur->next = prev;
        prev = cur;
        size = read(fd, cur, poolsize);
      } while( size == poolsize );
      if( size != (char *)&cur->region[state.ncur] - (char *)cur ) break;
    } StateReadClose(t, fd);
  }

  if( ini ) {
    Alloc(cur, 1);
    cur->next = NULL;
    state.ncur = t->nregions = 1;

    region = cur->region;
    region->div = 0;
    for( dim = 0; dim < t->ndim; ++dim ) {
      Bounds *b = &region->bounds[dim];
      b->lower = 0;
      b->upper = 1;
    }

    t->neval = 0;
    Sample(t, region);

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Totals *tot = &state.totals[comp];
      Result *r = &region->result[comp];
      tot->avg = tot->lastavg = tot->guess = r->avg;
      tot->err = tot->lasterr = r->err;
      tot->weightsum = 1/Max(Sq(r->err), NOTZERO);
      tot->avgsum = tot->weightsum*r->avg;
      tot->chisq = tot->chisqsum = tot->chisum = 0;
    }
  }

  /* main iteration loop */
  for( ; ; ) {
    count maxcomp, bisectdim;
    real maxratio, maxerr;
    Result result[NCOMP];
    Region *regionL, *regionR;
    Bounds *bL, *bR;

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        t->nregions, t->neval);

      for( comp = 0; comp < t->ncomp; ++comp ) {
        cTotals *tot = &state.totals[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, tot->avg, tot->err, tot->chisq, t->nregions - 1);
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
    regionL = cur->region;
    npool = state.ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region *region = &pool->region[ipool];
        creal err = region->result[maxcomp].err;
        if( err > maxerr ) {
          maxerr = err;
          regionL = region;
        }
      }

    if( state.ncur == POOLSIZE ) {
      Pool *prev = cur;
      Alloc(cur, 1);
      cur->next = prev;
      state.ncur = 0;
    }
    regionR = &cur->region[state.ncur++];

    regionR->div = ++regionL->div;
    FCopy(result, regionL->result);
    XCopy(regionR->bounds, regionL->bounds);

    bisectdim = result[maxcomp].bisectdim;
    bL = &regionL->bounds[bisectdim];
    bR = &regionR->bounds[bisectdim];
    bL->upper = bR->lower = .5*(bL->upper + bL->lower);

    Sample(t, regionL);
    Sample(t, regionR);

    for( comp = 0; comp < t->ncomp; ++comp ) {
      cResult *r = &result[comp];
      Result *rL = &regionL->result[comp];
      Result *rR = &regionR->result[comp];
      Totals *tot = &state.totals[comp];
      real diff, err, w, avg, sigsq;

      tot->lastavg += diff = rL->avg + rR->avg - r->avg;

      diff = fabs(.25*diff);
      err = rL->err + rR->err;
      if( err > 0 ) {
        creal c = 1 + 2*diff/err;
        rL->err *= c;
        rR->err *= c;
      }
      rL->err += diff;
      rR->err += diff;
      tot->lasterr += rL->err + rR->err - r->err;

      tot->weightsum += w = 1/Max(tot->lasterr, NOTZERO);
      sigsq = 1/tot->weightsum;
      tot->avgsum += w*tot->lastavg;
      avg = sigsq*tot->avgsum;
      tot->chisum += w *= tot->lastavg - tot->guess;
      tot->chisqsum += w*tot->lastavg;
      tot->chisq = tot->chisqsum - avg*tot->chisum;

      if( LAST ) {
        tot->avg = tot->lastavg;
        tot->err = tot->lasterr;
      }
      else {
        tot->avg = avg;
        tot->err = sqrt(sigsq);
      }
    }
    ++t->nregions;

    if( StateWriteTest(t) ) {
      StateWriteOpen(t, fd) {
        Pool *prev = cur;
        state.signature = StateSignature(t, 4);
        state.nregions = t->nregions;
        state.neval = t->neval;
        write(fd, &state, sizeof state);
        while( (prev = prev->next) ) write(fd, prev, poolsize);
fprintf(stderr, "ncur=%d  sizeof(Region)=%ld  size=%ld\n", state.ncur, 
sizeof(Region), (char *)&cur->region[state.ncur] - (char *)cur);
        write(fd, cur, (char *)&cur->region[state.ncur] - (char *)cur);
      } StateWriteClose(t, fd);
    }
  }

  for( comp = 0; comp < t->ncomp; ++comp ) {
    cTotals *tot = &state.totals[comp];
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, t->nregions - 1);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", t->nregions);

    npool = state.ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region const *region = &pool->region[ipool];
        real lower[NDIM], upper[NDIM];

        for( dim = 0; dim < t->ndim; ++dim ) {
          cBounds *b = &region->bounds[dim];
          lower[dim] = b->lower;
          upper[dim] = b->upper;
        }

        MLPutFunction(stdlink, "Cuba`Cuhre`region", 3);
        MLPutRealList(stdlink, lower, t->ndim);
        MLPutRealList(stdlink, upper, t->ndim);

        MLPutFunction(stdlink, "List", t->ncomp);
        for( comp = 0; comp < t->ncomp; ++comp ) {
          cResult *r = &region->result[comp];
          real res[] = {r->avg, r->err};
          MLPutRealList(stdlink, res, Elements(res));
        }
      }
  }
#endif

abort:
  while( (pool = cur) ) {
    cur = cur->next;
    free(pool);
  }

  WaitCores(t);
  FrameFree(t);
  RuleFree(t);

  StateRemove(t);

  return fail;
}

