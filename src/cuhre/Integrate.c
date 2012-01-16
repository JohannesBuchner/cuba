/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Cuhre
		last modified 8 Apr 09 th
*/


#define POOLSIZE 1024

static int Integrate(creal epsrel, creal epsabs,
  cint flags, number mineval, cnumber maxeval, ccount key,
  real *integral, real *error, real *prob)
{
  TYPEDEFREGION;
  typedef struct pool {
    struct pool *next;
    Region region[POOLSIZE];
  } Pool;

  count dim, comp, ncur, nregions, ipool, npool;
  int fail = 1;
  Rule rule;
  Totals totals[NCOMP];
  Pool *cur = NULL, *pool;
  Region *region;

  if( VERBOSE > 1 ) {
    char s[256];
    sprintf(s, "Cuhre input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key " COUNT,
      ndim_, ncomp_,
      epsrel, epsabs,
      flags, mineval, maxeval,
      key);
    Print(s);
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  if( key == 13 && ndim_ == 2 ) Rule13Alloc(&rule);
  else if( key == 11 && ndim_ == 3 ) Rule11Alloc(&rule);
  else if( key == 9 ) Rule9Alloc(&rule);
  else if( key == 7 ) Rule7Alloc(&rule);
  else {
    if( ndim_ == 2 ) Rule13Alloc(&rule);
    else if( ndim_ == 3 ) Rule11Alloc(&rule);
    else Rule9Alloc(&rule);
  }

  Alloc(rule.x, rule.n*(ndim_ + ncomp_));
  rule.f = rule.x + rule.n*ndim_;

  mineval = IMax(mineval, rule.n + 1);

  Alloc(cur, 1);
  cur->next = NULL;
  ncur = 1;

  region = cur->region;
  region->div = 0;
  for( dim = 0; dim < ndim_; ++dim ) {
    Bounds *b = &region->bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }

  Sample(&rule, region, flags);

  for( comp = 0; comp < ncomp_; ++comp ) {
    Totals *tot = &totals[comp];
    Result *r = &region->result[comp];
    tot->avg = tot->lastavg = tot->guess = r->avg;
    tot->err = tot->lasterr = r->err;
    tot->weightsum = 1/Max(Sq(r->err), NOTZERO);
    tot->avgsum = tot->weightsum*r->avg;
    tot->chisq = tot->chisqsum = tot->chisum = 0;
  }

  for( nregions = 1; ; ++nregions ) {
    count maxcomp, bisectdim;
    real maxratio, maxerr;
    Result result[NCOMP];
    Region *regionL, *regionR;
    Bounds *bL, *bR;

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        nregions, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cTotals *tot = &totals[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, tot->avg, tot->err, tot->chisq, nregions - 1);
      }

      Print(s);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( comp = 0; comp < ncomp_; ++comp ) {
      creal ratio = totals[comp].err/MaxErr(totals[comp].avg);
      if( ratio > maxratio ) {
        maxratio = ratio;
        maxcomp = comp;
      }
    }

    if( maxratio <= 1 && neval_ >= mineval ) {
      fail = 0;
      break;
    }

    if( neval_ >= maxeval ) break;

    maxerr = -INFTY;
    regionL = cur->region;
    npool = ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region *region = &pool->region[ipool];
        creal err = region->result[maxcomp].err;
        if( err > maxerr ) {
          maxerr = err;
          regionL = region;
        }
      }

    if( ncur == POOLSIZE ) {
      Pool *prev = cur;
      Alloc(cur, 1);
      cur->next = prev;
      ncur = 0;
    }
    regionR = &cur->region[ncur++];

    regionR->div = ++regionL->div;
    ResCopy(result, regionL->result);
    VecCopy(regionR->bounds, regionL->bounds);

    bisectdim = result[maxcomp].bisectdim;
    bL = &regionL->bounds[bisectdim];
    bR = &regionR->bounds[bisectdim];
    bL->upper = bR->lower = .5*(bL->upper + bL->lower);

    Sample(&rule, regionL, flags);
    Sample(&rule, regionR, flags);

    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &result[comp];
      Result *rL = &regionL->result[comp];
      Result *rR = &regionR->result[comp];
      Totals *tot = &totals[comp];
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

      tot->weightsum += w = 1/Max(Sq(tot->lasterr), NOTZERO);
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
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cTotals *tot = &totals[comp];
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, nregions - 1);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", nregions);

    npool = ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region const *region = &pool->region[ipool];
        real lower[NDIM], upper[NDIM];

        for( dim = 0; dim < ndim_; ++dim ) {
          cBounds *b = &region->bounds[dim];
          lower[dim] = b->lower;
          upper[dim] = b->upper;
        }

        MLPutFunction(stdlink, "Cuba`Cuhre`region", 3);
        MLPutRealList(stdlink, lower, ndim_);
        MLPutRealList(stdlink, upper, ndim_);

        MLPutFunction(stdlink, "List", ncomp_);
        for( comp = 0; comp < ncomp_; ++comp ) {
          cResult *r = &region->result[comp];
          real res[] = {r->avg, r->err};
          MLPutRealList(stdlink, res, Elements(res));
        }
      }
  }
#endif

#ifdef MLVERSION
abort:
#endif

  while( (pool = cur) ) {
    cur = cur->next;
    free(pool);
  }

  free(rule.x);
  RuleFree(&rule);

  nregions_ = nregions;

  return fail;
}

