#include "Sobol.c"
#include "ChiSquare.c"
#include "Grid.c"
#include "Sample.c"
#include "Fluct.c"
#include "Integrate.c"

#define MINDIM SOBOL_MINDIM

#define MAXDIM SOBOL_MAXDIM

#if NDIM > 0 && NDIM < MAXDIM
#undef MAXDIM
#define MAXDIM NDIM
#endif

