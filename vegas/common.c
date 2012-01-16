#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"
#include "Integrate.c"

static inline bool BadDimension(ccount ndim, cint flags)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
}

