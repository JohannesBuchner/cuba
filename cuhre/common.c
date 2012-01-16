#include "ChiSquare.c"
#include "Rule.c"
#include "Integrate.c"

static inline bool BadDimension(ccount ndim)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < 2;
}
