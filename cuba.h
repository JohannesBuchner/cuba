/*
	cuba.h
		Prototypes for the Cuba library
		this file is part of Cuba
		last modified 11 Dec 13 th
*/

#ifdef __cplusplus
extern "C" {
#endif

	/* integrand_t is intentionally a minimalistic integrand type.
	   It includes neither the nvec argument (for vectorized
	   evaluation) nor the extra arguments passed by Vegas/Suave
	   (weight, iter) and Divonne (phase).
	   In most cases, integrand_t is just what you want, otherwise
	   simply use an explicit typecast to integrand_t in the Cuba
	   invocation. */
typedef int (*integrand_t)(const int *ndim, const double x[],
  const int *ncomp, double f[], void *userdata);

typedef void (*peakfinder_t)(const int *ndim, const double b[],
  int *n, double x[]);

void Vegas(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int nstart, const int nincrease, const int nbatch,
  const int gridno, const char *statefile,
  int *neval, int *fail,
  double integral[], double error[], double prob[]);

void llVegas(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const long long int nstart, const long long int nincrease,
  const long long int nbatch,
  const int gridno, const char *statefile,
  long long int *neval, int *fail,
  double integral[], double error[], double prob[]);

void Suave(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int nnew, const double flatness,
  const char *statefile,
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);

void llSuave(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const long long int nnew, const double flatness,
  const char *statefile,
  int *nregions, long long int *neval, int *fail,
  double integral[], double error[], double prob[]);

void Divonne(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const int mineval, const int maxeval,
  const int key1, const int key2, const int key3, const int maxpass,
  const double border, const double maxchisq, const double mindeviation,
  const int ngiven, const int ldxgiven, double xgiven[],
  const int nextra, peakfinder_t peakfinder,
  const char *statefile,
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);

void llDivonne(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int seed,
  const long long int mineval, const long long int maxeval,
  const int key1, const int key2, const int key3, const int maxpass,
  const double border, const double maxchisq, const double mindeviation,
  const long long int ngiven, const int ldxgiven, double xgiven[],
  const long long int nextra,
  void (*peakfinder)(const int *, const double [], int *, double []),
  const char *statefile,
  int *nregions, long long int *neval, int *fail,
  double integral[], double error[], double prob[]);

void Cuhre(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const int nvec,
  const double epsrel, const double epsabs,
  const int flags, const int mineval, const int maxeval,
  const int key,
  const char *statefile,
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);

void llCuhre(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const double epsrel, const double epsabs,
  const int flags,
  const long long int mineval, const long long int maxeval,
  const int key,
  const char *statefile,
  int *nregions, long long int *neval, int *fail,
  double integral[], double error[], double prob[]);

void cubasetinit(void (*)(), void *);
void cubasetexit(void (*)(), void *);
void cubaruninit(void);
void cubaruninit(void);

#ifdef __cplusplus
}
#endif

