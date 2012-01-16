#define DEBFile stdout

#define DEBf "%9.5f"
//#define DEBf "%21.17f"
//#define DEBf "%a"


#ifdef DEBFile

#define DEBOpen()
#define DEBClose()

#else

FILE *DEBFile;

static inline void DEBOpen()
{
#ifdef MLVERSION
  DEBFile = fopen("log-mma", "w");
#else
  DEBFile = fopen("log-c", "w");
#endif
}

static inline void DEBClose()
{
  fclose(DEBFile);
}

#endif


#define DEB(...) fprintf(DEBFile, __VA_ARGS__); fflush(DEBFile)


static inline void DEBVec(const char *s, creal *d)
{
  char space[strlen(s) + 2];
  count dim;

  memset(space, ' ', sizeof(space));
  space[sizeof(space) - 1] = 0;

  DEB("%s=" DEBf "\n", s, d[0]);
  for( dim = 1; dim < ndim_; ++dim )
    DEB("%s" DEBf "\n", space, d[dim]);
}


/*
static inline void DEBRegion(const char *s, cBounds *b)
{
  char space[strlen(s) + 3];
  count dim;

  memset(space, ' ', sizeof(space));
  space[sizeof(space) - 1] = 0;

  DEB("%s: " DEBf " - " DEBf "\n", s, b[0].lower, b[0].upper);
  for( dim = 1; dim < ndim_; ++dim )
    DEB("%s" DEBf " - " DEBf "\n", space, b[dim].lower, b[dim].upper);
}
*/


static inline void DEBMem(const char *s)
{
  int kbytes = -1;
  FILE *f;
  char procfile[128];

  sprintf(procfile, "/proc/%d/status", getpid());
  f = fopen(procfile, "r");
  while( !feof(f) ) {
    char s[128];
    *s = 0;
    fgets(s, sizeof(s), f);
    if( sscanf(s, "VmSize: %d", &kbytes) == 1 ) break;
  }
  fclose(f);

  DEB("MEM %s: %d kbytes\n", s, kbytes);
}

