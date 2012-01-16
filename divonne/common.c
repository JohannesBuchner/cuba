static void Explore(void *voidregion, cSamples *samples, cint depth, cint flags);

static void Split(void *voidregion, int depth);

#include "Sobol.c"
#include "ChiSquare.c"
#include "Rule.c"
#include "Sample.c"
#include "FindMinimum.c"
#include "Explore.c"
#include "Split.c"
#include "Integrate.c"

#if KOROBOV_MINDIM > SOBOL_MINDIM
#define MINDIM KOROBOV_MINDIM
#else
#define MINDIM SOBOL_MINDIM
#endif

#if KOROBOV_MAXDIM < SOBOL_MAXDIM
#define MAXDIM KOROBOV_MAXDIM
#else
#define MAXDIM SOBOL_MAXDIM
#endif

#if NDIM > 0 && NDIM < MAXDIM
#undef MAXDIM
#define MAXDIM NDIM
#endif

