/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 19 Dec 11 th
*/


#include "ChiSquare.c"
#include "Rule.c"

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > NDIM ) return true;
  return t->ndim < 2;
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

