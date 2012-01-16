/*
	DoSample.c
		the actual sampling routine, serial and parallel,
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 24 Nov 11 th
*/

#define MINSLICE 10
#define MINCORES 1
//#define MINCORES 2

#if defined(VEGAS) || defined(SUAVE)
#define VEG_ONLY(...) __VA_ARGS__
#else
#define VEG_ONLY(...)
#endif

#ifdef DIVONNE
#define DIV_ONLY(...) __VA_ARGS__
#define LDX(ldx) ldx
#else
#define DIV_ONLY(...)
#define LDX(ldx) t->ndim
#endif

typedef struct {
  real *f;
  number n;
  VEG_ONLY(count iter;)
  DIV_ONLY(number neval_opt, neval_cut;
           count ldx, phase, iregion;)
#define NREGIONS ldx
#define NEVAL n
#define RETVAL phase
} Slice;

/*********************************************************************/

#ifndef MSG_WAITALL
/* Windows */
#define MSG_WAITALL 0
#endif

static inline int readsock(int fd, void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = recv(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

static inline int writesock(int fd, const void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = send(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

/*********************************************************************/

static inline bool SampleSerial(cThis *t, number n, creal *x, real *f
  VEG_ONLY(, creal *w, ccount iter)
  DIV_ONLY(, ccount ldx))
{
  while( n-- ) {
    if( t->integrand(&t->ndim, x, &t->ncomp, f, t->userdata
          VEG_ONLY(, w++, &iter)
          DIV_ONLY(, &t->phase)) == ABORT ) return true;
    x += LDX(ldx);
    f += t->ncomp;
  }
  return false;
}

/*********************************************************************/

static void DoSample(This *t, number n, creal *x, real *f
  VEG_ONLY(, creal *w, ccount iter)
  DIV_ONLY(, ccount ldx))
{
  char s[128];
  Slice slice;
  int ncores;

  t->neval += n;

  ncores = IMin(t->ncores, n/MINSLICE);

  if( ncores < MINCORES ) {
    if( VERBOSE > 2 ) {
      sprintf(s, "sampling " NUMBER " points serially", n);
      Print(s);
    }

    if( SampleSerial(t, n, x, f
          VEG_ONLY(, w, iter)
          DIV_ONLY(, ldx)) ) longjmp(t->abort, -99);
  }
  else {
    int core, abort;

    slice.n = (n + ncores - 1)/ncores;

    if( VERBOSE > 2 ) {
      sprintf(s, "sampling " NUMBER " points each on %d cores",
        slice.n, ncores);
      Print(s);
    }

    slice.f = f;
    VEG_ONLY(slice.iter = iter;)
    DIV_ONLY(slice.ldx = ldx;)
    DIV_ONLY(slice.phase = t->phase;)

    for( core = 0; core < ncores; ++core ) {
      cint fd = t->child[core];
      writesock(fd, &slice, sizeof slice);
      VEG_ONLY(writesock(fd, w, slice.n*sizeof *w);)
      writesock(fd, x, slice.n*LDX(ldx)*sizeof *x);

      VEG_ONLY(w += n;)
      x += slice.n*LDX(ldx);
      slice.f += slice.n*t->ncomp;
      n -= slice.n;
      slice.n = IMin(slice.n, n);
    }

    abort = 0;
    for( core = ncores; --core >= 0; ) {
      cint fd = t->child[core];
      readsock(fd, &slice, sizeof slice);
      if( slice.n == 0 ) abort = 1;
      else readsock(fd, slice.f, slice.n*t->ncomp*sizeof *f);
    }
    if( abort ) longjmp(t->abort, -99);
  }
}

/*********************************************************************/

#ifdef DIVONNE
static inline int ReadyCore(cThis *t)
{
  int core;
  fd_set ready;

  memcpy(&ready, &t->children, sizeof ready);
  select(t->nchildren, &ready, NULL, NULL, NULL);

  for( core = 0; core < t->ncores; ++core )
    if( FD_ISSET(t->child[core], &ready) ) break;

  return core;
}

/*********************************************************************/

static int ExploreParent(This *t, cint iregion)
{
  TYPEDEFREGION;
  Region *region;
  Slice slice;
  int ireg = iregion, core = t->running;

  if( t->ncores < MINCORES ) return Explore(t, iregion);

  if( t->running >= ((iregion < 0) ? 1 : t->ncores) ) {
    Totals totals[t->ncomp];
    count comp, succ;
    cint fd = t->child[core = ReadyCore(t)];

    --t->running;
    readsock(fd, &slice, sizeof slice);
//DEBSLICE("parent read", fd, slice);
    ireg = slice.iregion;
    region = RegionPtr(ireg);
    succ = ireg + region->next;
    readsock(fd, region, sizeof(Region));
    if( --slice.NREGIONS > 0 ) {
      region->next = t->nregions - ireg;
      EnlargeRegions(t, slice.NREGIONS);
      readsock(fd, RegionPtr(t->nregions), slice.NREGIONS*sizeof(Region));
      t->nregions += slice.NREGIONS;
      RegionPtr(t->nregions-1)->next = succ - t->nregions + 1;
    }

    readsock(fd, totals, sizeof totals);
    for( comp = 0; comp < t->ncomp; ++comp )
      t->totals[comp].secondspread =
        Max(t->totals[comp].secondspread, totals[comp].secondspread);

    t->neval += slice.NEVAL;
    t->neval_opt += slice.neval_opt;
    t->neval_cut += slice.neval_cut;

    if( slice.RETVAL == -1 ) return -1;
  }

  if( iregion >= 0 ) {
    region = RegionPtr(iregion);
    cint fd = t->child[core];
    slice.n = 0;
    slice.phase = t->phase;
    slice.iregion = iregion;
//DEBSLICE("  parent write", fd, slice);
    writesock(fd, &slice, sizeof slice);
    writesock(fd, &t->samples[region->isamples], sizeof(Samples));
    writesock(fd, region, sizeof *region);
    writesock(fd, t->totals, sizeof *t->totals);
    region->depth = 0;
    ++t->running;
  }

  return ireg;
}
#endif

/*********************************************************************/

static inline void DoChild(This *t, cint fd)
{
  Slice slice;

#ifdef DIVONNE
  TYPEDEFREGION;
  Totals totals[t->ncomp];

  t->totals = totals;
  t->ncores = 0;	/* no recursive forks */
  AllocRegions(t);
  SamplesIni(&t->samples[0]);
  t->samples[0].n = 0;
  SamplesIni(&t->samples[1]);
  t->samples[1].n = 0;
  SamplesIni(&t->samples[2]);
  t->samples[2].n = 0;
#endif

  while( readsock(fd, &slice, sizeof slice) ) {
    number n = slice.n;
    DIV_ONLY(t->phase = slice.phase;)
//DEBSLICE("  child read", fd, slice);
    if( n > 0 ) {
      VEG_ONLY(real w[n];)
      real x[n*LDX(slice.ldx)];
      real f[n*t->ncomp];

      VEG_ONLY(readsock(fd, w, sizeof w);)
      readsock(fd, x, sizeof x);

      if( SampleSerial(t, n, x, f
            VEG_ONLY(, w, slice.iter)
            DIV_ONLY(, slice.ldx)) ) slice.n = 0;
      writesock(fd, &slice, sizeof slice);
      if( slice.n ) writesock(fd, f, sizeof f);
    }
#ifdef DIVONNE
    else {
      Samples *samples, psamples;

      readsock(fd, &psamples, sizeof psamples);
      readsock(fd, RegionPtr(0), sizeof(Region));
      readsock(fd, totals, sizeof totals);
      t->nregions = 1;
      t->neval = t->neval_opt = t->neval_cut = 0;

      samples = &t->samples[RegionPtr(0)->isamples];
      if( psamples.n != samples->n ) {
        SamplesFree(samples);
        *samples = psamples;
        SamplesAlloc(t, samples);
      }

      slice.RETVAL = Explore(t, 0);
      slice.NREGIONS = t->nregions;
      slice.NEVAL = t->neval;
      slice.neval_opt = t->neval_opt;
      slice.neval_cut = t->neval_cut;
//DEBSLICE("child write", fd, slice);
      writesock(fd, &slice, sizeof slice);
      writesock(fd, RegionPtr(0), t->nregions*sizeof(Region));
      writesock(fd, totals, sizeof totals);
    }
#endif
  }

  exit(0);
}

/*********************************************************************/

static inline void ForkCores(This *t)
{
  int core;
  cchar *env = getenv("CUBACORES");

  t->ncores = env ? atoi(env) : sysconf(_SC_NPROCESSORS_ONLN);
#ifdef HAVE_GETLOADAVG
  if( env == NULL || t->ncores < 0 ) {
    double load = 0;
    getloadavg(&load, 1);
    t->ncores = abs(t->ncores) - floor(load);
  }
#endif

#ifdef DIVONNE
  t->nchildren = t->running = 0;
#endif

  if( t->ncores < MINCORES ) return;
  if( VERBOSE ) printf("using %d cores\n", t->ncores);
  fflush(stdout);

  Alloc(t->child, t->ncores);
  for( core = 0; core < t->ncores; ++core ) {
    int fd[2];
    pid_t pid;
    assert(
      socketpair(AF_LOCAL, SOCK_STREAM, 0, fd) != -1 &&
      (pid = fork()) != -1 );
    if( pid == 0 ) {
      close(fd[0]);
      DoChild(t, fd[1]);
    }
    close(fd[1]);
    t->child[core] = fd[0];
#ifdef DIVONNE
    FD_SET(fd[0], &t->children);
    t->nchildren = IMax(t->nchildren, fd[0] + 1);
#endif
  }
}

/*********************************************************************/

static inline void WaitCores(cThis *t)
{
  if( t->ncores >= MINCORES ) {
    int core;
    pid_t pid;
    for( core = 0; core < t->ncores; ++core )
      close(t->child[core]);
    free(t->child);
    for( core = 0; core < t->ncores; ++core )
      wait(&pid);
  }
}

