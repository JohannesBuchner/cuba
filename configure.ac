# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([Cuba], [1.4], [hahn@feynarts.de])
AC_CONFIG_SRCDIR([cuba.h])

LIBS=-lm
USER_CFLAGS="$CFLAGS"

AC_PROG_CC
AC_PROG_F77
AC_PROG_RANLIB
AC_PROG_INSTALL

if test $GCC = yes ; then
  CFLAGS=${USER_CFLAGS:--O3 -fomit-frame-pointer -ffast-math}
else
  CFLAGS=${USER_CFLAGS:--O}
fi

AC_C_CONST
AC_C_INLINE
AC_C_LONG_DOUBLE

AC_CHECK_FUNCS([powl])
AC_CHECK_FUNCS([erf])

MAXDIM=${MAXDIM:-16}
AC_ARG_WITH(maxdim,
  [AS_HELP_STRING([--with-maxdim=N],
    [the maximum dimension for integration,
     if variable-size array are not supported])],
  [MAXDIM=$withval])

MAXCOMP=${MAXCOMP:-4}
AC_ARG_WITH(maxcomp,
  [AS_HELP_STRING([--with-maxcomp=N],
    [the maximum number of components of the integrand,
     if variable-size array are not supported])],
  [MAXCOMP=$withval])

AC_MSG_CHECKING([for variable-size arrays])
AC_COMPILE_IFELSE([AC_LANG_SOURCE(,[[
	void test(int n)
	{
	  char s[n];
	}
  ]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no, using MAXDIM=$MAXDIM and MAXCOMP=$MAXCOMP])
   AC_DEFINE_UNQUOTED([NDIM], [$MAXDIM], [Maximum number of components])
   AC_DEFINE_UNQUOTED([NCOMP], [$MAXCOMP], [Maximum number of dimensions])]
)

AC_CHECK_PROGS([HAVE_QMAKE], [qmake])
test -n "$HAVE_QMAKE" && TOOLS_DEFAULT=tools
AC_SUBST([TOOLS_DEFAULT])

test -n "$F77" && DEMO_FORTRAN_DEFAULT=demo-fortran
AC_SUBST([DEMO_FORTRAN_DEFAULT])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([makefile])

AC_OUTPUT

