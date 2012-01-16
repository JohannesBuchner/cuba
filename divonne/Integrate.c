/*
	Integrate.c
		partition the integration region until each region
		has approximately equal spread = 1/2 vol (max - min),
		then do a final integration over all regions
		this file is part of Divonne
		last modified 21 Jan 05 th
*/


#define INIDEPTH 3
#define DEPTH 5

/*********************************************************************/

static int Integrate(creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  int key1, int key2, int key3, ccount maxpass, 
  creal maxchisq, creal mindeviation,
  real *integral, real *error, real *prob)
{
  TYPEDEFREGION;

  Region anchor, *region;
  Totals totals[NCOMP];
  real nneed, weight;
  count dim, comp, iter, pass = 0, err;
  number nwant, nmin = INT_MAX;
  int fail = -1;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Divonne input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key1 %d\n  key2 %d\n  key3 %d\n  maxpass " COUNT "\n"
      "  border " REAL "\n  maxchisq " REAL "\n  mindeviation " REAL "\n"
      "  ngiven " NUMBER "\n  nextra " NUMBER "\n",
      ndim_, ncomp_,
      epsrel, epsabs,
      flags, mineval, maxeval,
      key1, key2, key3, maxpass,
      border_.lower, maxchisq, mindeviation,
      ngiven_, nextra_);
    Print(s);
  }

  anchor.next = NULL;
  for( dim = 0; dim < ndim_; ++dim ) {
    Bounds *b = &anchor.bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }

  RuleIni(&rule7_);
  RuleIni(&rule9_);
  RuleIni(&rule11_);
  RuleIni(&rule13_);
  SamplesIni(&samples_[0]);
  SamplesIni(&samples_[1]);
  SamplesIni(&samples_[2]);

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  /* Step 1: partition the integration region */

  if( VERBOSE ) Print("Partitioning phase:");

  if( IsSobol(key1) || IsSobol(key2) || IsSobol(key3) )
    IniRandom(2*maxeval, flags);

  SamplesLookup(&samples_[0], key1,
    (number)47, (number)INT_MAX, (number)0);
  SamplesAlloc(&samples_[0]);

  totals_ = totals;
  Zap(totals);
  phase_ = 1;

  Explore(&anchor, &samples_[0], INIDEPTH, 1);

  for( iter = 1; ; ++iter ) {
    Totals *maxtot;

    for( comp = 0; comp < ncomp_; ++comp ) {
      Totals *tot = &totals[comp];
      tot->avg = tot->spread = 0;
      tot->maxspread = tot->secondspread = -INFTY;
    }

    nregions_ = 0;
    for( region = anchor.next; region; region = region->next ) {
      ++nregions_;
      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];
        tot->avg += r->avg;
        tot->spread += Sq(r->spread);
        if( r->spread > tot->maxspread ) {
          tot->secondspread = tot->maxspread;
          tot->maxspread = r->spread;
          tot->maxregion = region;
        }
        else if( r->spread > tot->secondspread )
          tot->secondspread = r->spread;
      }
    }

    maxtot = totals;
    for( comp = 0; comp < ncomp_; ++comp ) {
      Totals *tot = &totals[comp];
      integral[comp] = tot->avg;
      tot->spread = sqrt(tot->spread);
      if( tot->spread > maxtot->spread ) maxtot = tot;
      error[comp] = tot->spread*samples_[0].weight;
    }

    if( VERBOSE ) {
      char s[128 + 64*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " COUNT " regions\n"
        NUMBER7 " integrand evaluations so far,\n"
        NUMBER7 " in optimizing regions,\n"
        NUMBER7 " in finding cuts",
        iter, nregions_, neval_, neval_opt_, neval_cut_);

      for( comp = 0; comp < ncomp_; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL,
          comp + 1, integral[comp], error[comp]);

      Print(s);
    }

    if( neval_ > maxeval ) break;

    nneed = maxtot->spread/Max(fabs(maxtot->avg)*epsrel, epsabs);
    if( nneed < MAXPRIME ) {
      cnumber n = neval_ + nregions_*(number)(nneed + .5);
      if( n < nmin ) {
        nmin = n;
        pass = 0;
      }
      else if( ++pass > maxpass && n >= mineval ) break;
    }

    Split(maxtot->maxregion, DEPTH);
  }

  /* Step 2: do a "full" integration on each region */

/* nneed = samples_[0].neff + 1; */
  nneed = 2*samples_[0].neff;
  for( comp = 0; comp < ncomp_; ++comp ) {
    Totals *tot = &totals[comp];
    creal maxspread = Max(fabs(tot->avg)*epsrel, epsabs);
    nneed = Max(nneed, tot->spread /= maxspread);
    tot->maxspread = Sq(mindeviation*maxspread);
  }
  nwant = Min(nneed + .5, MARKMASK/40);

  err = SamplesLookup(&samples_[1], key2, nwant,
    (maxeval - neval_)/nregions_ + 1, samples_[0].n + 1);

  /* the number of points needed to reach the desired accuracy */
  fail = Unmark(err)*nregions_;

  if( Marked(err) ) {
    if( VERBOSE ) Print("\nNot enough samples left for final integration.");
    for( comp = 0; comp < ncomp_; ++comp )
      prob[comp] = -999;
    weight = samples_[0].weight;
  }
  else {
    count df, nlimit;

    SamplesAlloc(&samples_[1]);

    if( VERBOSE ) {
      char s[128];
      sprintf(s, "\nFinal integration on " COUNT " regions with " NUMBER " samples per region.",
        nregions_, samples_[1].neff);
      Print(s);
    }

    ResClear(integral);
    ResClear(error);
    ResClear(prob);

    nlimit = maxeval - nregions_*samples_[1].n;
    df = nregions_ = 0;

    for( region = anchor.next; region; region = region->next ) {
      char s[64*NDIM + 256*NCOMP], *p = s;
      int todo;

refine:
      phase_ = 2;
      samples_[1].sampler(&samples_[1], region->bounds, region->vol);
      nlimit += samples_[1].n;
      todo = 0;

      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];

        samples_[0].avg[comp] = r->avg;
        samples_[0].err[comp] = r->err;

        if( neval_ < nlimit ) {
          creal avg2 = samples_[1].avg[comp];
          creal err2 = samples_[1].err[comp];
          creal diffsq = Sq(avg2 - r->avg);

#define Var(s) Sq((s.err[comp] == 0) ? r->spread*s.weight : s.err[comp])

          if( err2*tot->spread > r->spread ||
              diffsq > Max(maxchisq*(Var(samples_[0]) + Var(samples_[1])),
                           EPS*Sq(avg2)) ) {
            if( key3 && diffsq > tot->maxspread ) {
              if( key3 == 1 ) {
                if( VERBOSE > 2 ) Print("\nSplit");
                phase_ = 1;
                Explore(region, &samples_[1], 1, 2);
                goto refine;
              }
              todo |= 3;
            }
            todo |= 1;
          }
        }
      }

      switch( todo ) {
      case 1:	/* get spread right */
        Explore(region, &samples_[1], 0, 2);
        break;

      case 3:	/* sample region again with more points */
        if( MEM(&samples_[2]) == NULL ) {
          SamplesLookup(&samples_[2], key3,
            nwant, (number)INT_MAX, (number)0);
          SamplesAlloc(&samples_[2]);
        }
        phase_ = 3;
        samples_[2].sampler(&samples_[2], region->bounds, region->vol);
        Explore(region, &samples_[2], 0, 2);
        ++region->depth;	/* misused for df here */
        ++df;
      }

      ++region->depth;	/* misused for df here */
      ++nregions_;

      if( VERBOSE > 2 ) {
        for( dim = 0; dim < ndim_; ++dim ) {
          cBounds *b = &region->bounds[dim];
          p += sprintf(p,
            (dim == 0) ? "\nRegion (" REALF ") - (" REALF ")" :
                         "\n       (" REALF ") - (" REALF ")",
            b->lower, b->upper);
        }
      }

      for( comp = 0; comp < ncomp_; ++comp ) {
        Result *r = &region->result[comp];

        creal x1 = samples_[0].avg[comp];
        creal s1 = Var(samples_[0]);
        creal x2 = samples_[1].avg[comp];
        creal s2 = Var(samples_[1]);
        creal r2 = (s1 == 0) ? Sq(samples_[1].neff*samples_[0].weight) : s2/s1;

        real norm = 1 + r2;
        real avg = x2 + r2*x1;
        real sigsq = s2;
        real chisq = Sq(x2 - x1);
        real chiden = s1 + s2;

        if( todo == 3 ) {
          creal x3 = samples_[2].avg[comp];
          creal s3 = Var(samples_[2]);
          creal r3 = (s2 == 0) ? Sq(samples_[2].neff*samples_[1].weight) : s3/s2;

          norm = 1 + r3*norm;
          avg = x3 + r3*avg;
          sigsq = s3;
          chisq = s1*Sq(x3 - x2) + s2*Sq(x3 - x1) + s3*chisq;
          chiden = s1*s2 + s3*chiden;
        }

        avg = LAST ? r->avg : (sigsq *= norm = 1/norm, avg*norm);
        if( chisq > EPS ) chisq /= Max(chiden, NOTZERO);

#define Out(s) s.avg[comp], r->spread*s.weight, s.err[comp]

        if( VERBOSE > 2 ) {
          p += sprintf(p, "\n[" COUNT "] "
            REAL " +- " REAL "(" REAL ")\n    "
            REAL " +- " REAL "(" REAL ")",
            comp + 1, Out(samples_[0]), Out(samples_[1]));
          if( todo == 3 ) p += sprintf(p, "\n    "
            REAL " +- " REAL "(" REAL ")",
            Out(samples_[2]));
          p += sprintf(p, "  \tchisq " REAL, chisq);
        }

        integral[comp] += avg;
        error[comp] += sigsq;
        prob[comp] += chisq;

        r->avg = avg;
        r->spread = sqrt(sigsq);
        r->chisq = chisq;
      }

      if( VERBOSE > 2 ) Print(s);
    }

    for( comp = 0; comp < ncomp_; ++comp )
      error[comp] = sqrt(error[comp]);

    df += nregions_;

    if( VERBOSE > 2 ) {
      char s[16 + 128*NCOMP], *p = s;

      p += sprintf(p, "\nTotals:");

      for( comp = 0; comp < ncomp_; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, integral[comp], error[comp], prob[comp], df);

      Print(s);
    }

    for( comp = 0; comp < ncomp_; ++comp )
      prob[comp] = ChiSquare(prob[comp], df);

    weight = 1;
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", nregions_);
    for( region = anchor.next; region; region = region->next ) {
      cBounds *b = region->bounds;
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < ndim_; ++dim ) {
        lower[dim] = b[dim].lower;
        upper[dim] = b[dim].upper;
      }

      MLPutFunction(stdlink, "Divonne`Private`region", 4);

      MLPutRealList(stdlink, lower, ndim_);
      MLPutRealList(stdlink, upper, ndim_);

      MLPutFunction(stdlink, "List", ncomp_);
      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->spread*weight, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->depth);  /* misused for df */
    }
  }
#endif

abort:
  SamplesFree(&samples_[2]);
  SamplesFree(&samples_[1]);
  SamplesFree(&samples_[0]);
  RuleFree(&rule13_);
  RuleFree(&rule11_);
  RuleFree(&rule9_);
  RuleFree(&rule7_);

  for( region = anchor.next; region; ) {
    Region *next = region->next;
    free(region);
    region = next;
  }

  return fail;
}

