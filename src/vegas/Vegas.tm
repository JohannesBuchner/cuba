:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Vegas::usage = "Vegas[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: NStart::usage = "NStart is an option of Vegas.
	It specifies the number of integrand evaluations per iteration to start with."

:Evaluate: NIncrease::usage = "NIncrease is an option of Vegas.
	It specifies the increase in the number of integrand evaluations per iteration."

:Evaluate: NBatch::usage = "NBatch is an option of Vegas.
	It specifies how many points are sent in one MathLink packet to be sampled by Mathematica."

:Evaluate: GridNo::usage = "GridNo is an option of Vegas.
	Vegas maintains an internal table in which it can memorize up to 10 grids, to be used on subsequent integrations.
	A GridNo between 1 and 10 selects the slot in this internal table.
	For other values the grid is initialized from scratch and discarded at the end of the integration."

:Evaluate: StateFile::usage = "StateFile is an option of Vegas.
	It specifies a file in which the internal state is stored after each iteration and from which it can be restored on a subsequent run.
	The state file is removed once the prescribed accuracy has been reached."

:Evaluate: MinPoints::usage = "MinPoints is an option of Vegas.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Vegas.
	It can take the values Last or All which determine whether only the last (largest) or all of the samples collected on a subregion over the iterations contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Vegas.
	If set to True, pseudo-random numbers are used instead of Sobol quasi-random numbers.
	If set to an integer value, that value is used as the seed for the pseudo-random number generator."

:Evaluate: SharpEdges::usage = "SharpEdges is an option of Vegas.
	It turns off smoothing of the importance function for integrands with sharp edges."

:Evaluate: $Weight::usage = "$Weight is a global variable set by Vegas during the evaluation of the integrand to the weight of the point being sampled."

:Evaluate: MapSample::usage = "MapSample is a function used to map the integrand over the points to be sampled."


:Evaluate: Begin["`Vegas`"]

:Begin:
:Function: Vegas
:Pattern: MLVegas[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, mineval_, maxeval_,
  nstart_, nincrease_,
  seed_, nbatch_, gridno_, state_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, mineval, maxeval,
  nstart, nincrease,
  seed, nbatch, gridno, state}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer, Integer,
  Integer, Integer,
  Integer, Integer, Integer, String}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Vegas] = {HoldFirst}

:Evaluate: Options[Vegas] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000,
	NStart -> 1000, NIncrease -> 500,
	NBatch -> 1000, GridNo -> 0, StateFile -> "",
	Verbose -> 1, Final -> All, PseudoRandom -> False,
	SharpEdges -> False, Compiled -> True}

:Evaluate: Vegas[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}],
	tags, vars, lower, range, integrand,
	rel, abs, mineval, maxeval, nstart, nincrease, nbatch,
	gridno, verbose, final, pseudo, edges, compiled},
	  Message[Vegas::optx, #, Vegas]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Vegas]];
	  {rel, abs, mineval, maxeval, nstart, nincrease, nbatch,
	    gridno, state, verbose, final, pseudo, edges, compiled} =
	    tags /. {opt} /. Options[Vegas];
	  {vars, lower, range} = Transpose[{v}];
	  range -= lower;
	  define[compiled, vars, lower, range, Simplify[Times@@ range]];
	  integrand = fun[f];
	  MLVegas[ndim, ncomp[f], 10.^-rel, 10.^-abs,
	    Min[Max[verbose, 0], 3] +
	      If[final === Last, 4, 0] +
	      If[pseudo =!= False, 8, 0] +
	      If[TrueQ[edges], 16, 0],
	    mineval, maxeval, nstart, nincrease,
	    If[IntegerQ[pseudo], pseudo, 0], nbatch, gridno, state]
	]

:Evaluate: Attributes[ncomp] = Attributes[fun] = {HoldAll}

:Evaluate: ncomp[f_List] := Length[f]

:Evaluate: _ncomp = 1

:Evaluate: define[True, vars:{v___}, lower_, range_, jac_] :=
	fun[f_] := Compile[{{t, _Real, 1}, w},
	  Block[{v, $Weight = w},
	    vars = lower + range t;
	    check[vars, Chop[f jac]//N] ]]

:Evaluate: define[_, vars:{v___}, lower_, range_, jac_] :=
	fun[f_] := Function[{t, w},
	  Block[{v, $Weight = w},
	    vars = lower + range t;
	    check[vars, Chop[f jac]//N] ]]

:Evaluate: check[_, f_Real] = {f}

:Evaluate: check[_, f:{__Real}] = f

:Evaluate: check[x_, _] := (Message[Vegas::badsample, ff, x]; {})

:Evaluate: sample[x_, w_] := Check[ Flatten @
	MapSample[integrand@@ # &, Transpose[{Partition[x, ndim], w}]],
	{} ]

:Evaluate: MapSample = Map

:Evaluate: Vegas::badsample = "`` is not a real-valued function at ``."

:Evaluate: Vegas::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Vegas::badcomp = "Cannot integrate `` components."

:Evaluate: Vegas::accuracy =
	"Desired accuracy was not reached within `` function evaluations."

:Evaluate: Vegas::success = "Needed `` function evaluations."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Vegas.tm
		Vegas Monte Carlo integration
		by Thomas Hahn
		last modified 5 Dec 08 th
*/


#include <setjmp.h>
#include "mathlink.h"
#include "util.c"

jmp_buf abort_;

/*********************************************************************/

static void Status(MLCONST char *msg, cint n)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 2);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Vegas");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n);
}

/*********************************************************************/

static void Print(MLCONST char *s)
{
  int pkt;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Print", 1);
  MLPutString(stdlink, s);
  MLEndPacket(stdlink);

  do {
    pkt = MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  } while( pkt != RETURNPKT );
}

/*********************************************************************/

static void DoSample(cnumber n, real *w, real *x, real *f)
{
  int pkt;
  real *mma_f;
  long mma_n;

  if( MLAbort ) goto abort;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`Vegas`sample", 2);
  MLPutRealList(stdlink, x, n*ndim_);
  MLPutRealList(stdlink, w, n);
  MLEndPacket(stdlink);

  while( (pkt = MLNextPacket(stdlink)) && (pkt != RETURNPKT) )
    MLNewPacket(stdlink);

  if( !MLGetRealList(stdlink, &mma_f, &mma_n) ) {
    MLClearError(stdlink);
    MLNewPacket(stdlink);
abort:
    MLPutFunction(stdlink, "Abort", 0);
    longjmp(abort_, 1);
  }

  if( mma_n != n*ncomp_ ) {
    MLDisownRealList(stdlink, mma_f, mma_n);
    MLPutSymbol(stdlink, "$Failed");
    longjmp(abort_, 1);
  }

  Copy(f, mma_f, n*ncomp_);
  MLDisownRealList(stdlink, mma_f, mma_n);

  neval_ += n;
}

/*********************************************************************/

#include "common.c"

void Vegas(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  cint seed, cint nbatch, cint gridno, const char *state)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) ) {
    Status("badcomp", ncomp);
    MLPutSymbol(stdlink, "$Failed");
  }
  else if( BadDimension(ndim, flags) ) {
    Status("baddim", ndim);
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    real integral[NCOMP], error[NCOMP], prob[NCOMP];
    count comp;
    int fail;

    neval_ = 0;
    mersenneseed = seed;
    vegasnbatch = nbatch;
    vegasgridno = gridno;
    strncpy(vegasstate, state, sizeof(vegasstate) - 1);
    vegasstate[sizeof(vegasstate) - 1] = 0;

    fail = Integrate(epsrel, epsabs,
      flags, mineval, maxeval, nstart, nincrease,
      integral, error, prob);

    Status(fail ? "accuracy" : "success", neval_);

    MLPutFunction(stdlink, "List", ncomp);
    for( comp = 0; comp < ncomp; ++comp ) {
      real res[] = {integral[comp], error[comp], prob[comp]};
      MLPutRealList(stdlink, res, Elements(res));
    }
  }

  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

