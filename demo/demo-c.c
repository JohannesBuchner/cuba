#include <stdio.h>
#include <math.h>
#include "cuba.h"


static inline double Sq(double x)
{
  return x*x;
}


static void Integrand(const int *ndim, const double xx[],
                      const int *ncomp, double ff[])
{
#define x xx[0]
#define y xx[1]
#define z xx[2]
#define f ff[0]

#define FUN 1

  const double rsq = Sq(x) + Sq(y) + Sq(z);

#if FUN == 1
  f = sin(x)*cos(y)*exp(z);
#elif FUN == 2
  f = 1/(Sq(x + y) + .003)*cos(y)*exp(z);
#elif FUN == 3
  f = 1/(3.75 - cos(M_PI*x) - cos(M_PI*y) - cos(M_PI*z));
#elif FUN == 4
  f = fabs(rsq - .125);
#elif FUN == 5
  f = exp(-rsq);
#elif FUN == 6
  f = 1/(1 - x*y*z + 1e-10);
#elif FUN == 7
  f = sqrt(fabs(x - y - z));
#elif FUN == 8
  f = exp(-x*y*z);
#elif FUN == 9
  f = Sq(x)/(cos(x + y + z + 1) + 5);
#elif FUN == 10
  f = (x > .5) ? 1/sqrt(x*y*z + 1e-5) : sqrt(x*y*z);
#else
  f = (rsq < 1) ? 1 : 0;
#endif
}

/*********************************************************************/

#define NDIM 3
#define NCOMP 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000

#define NSTART 1000
#define NINCREASE 500

#define NNEW 1000
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main()
{
  int comp, nregions, neval, fail;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];

  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand,
    EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
    NSTART, NINCREASE,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral[comp], error[comp], prob[comp]);

  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand,
    EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
    NNEW, FLATNESS,
    &nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral[comp], error[comp], prob[comp]);

  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand,
    EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
    KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    &nregions, &neval, &fail, integral, error, prob);

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral[comp], error[comp], prob[comp]);

  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand,
    EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
    KEY,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral[comp], error[comp], prob[comp]);

  return 0;
}

