CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -ffast-math

# uncomment the following if experiencing problems with a non-GNU cc,
# or if the code segfaults
#CFLAGS += -DNDIM=8 -DNCOMP=2

# if long double and/or powl are not available (even though C99
# requires them), uncomment the following line
#CFLAGS += -DNO_LONG_DOUBLE

HEADERS = cuba.h

LIB = libcuba.a

DEMO = demo-c demo-fortran

MMA = Vegas Suave Divonne Cuhre

all: lib demo mma

lib: $(LIB)

demo: $(DEMO)

mma: $(MMA)


PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

install:
	-mkdir -p $(LIBDIR) $(INCLUDEDIR) $(BINDIR)
	install -m 644 $(LIB) $(LIBDIR)
	install -m 644 $(HEADERS) $(INCLUDEDIR)
	install -s $(MMA) $(BINDIR)


check: demo
	./demo-c | grep RESULT > demo-c.out
	diff demo-c.out demo-c.out.ok


VEGAS_C = vegas/Vegas.c
VEGAS_F1 = vegas/vegas-f.c
VEGAS_F2 = vegas/vegas_f.c
VEGAS_MMA = vegas/Vegas.tm
VEGAS_H = vegas/decl.h vegas/stddecl.h
VEGAS_SRCS = $(VEGAS_H) vegas/util.c vegas/debug.c vegas/common.c \
  vegas/Random.c vegas/ChiSquare.c vegas/Grid.c vegas/Integrate.c

$(LIB)(Vegas.o): $(VEGAS_C) $(VEGAS_SRCS)
	$(CC) $(CFLAGS) -c -o Vegas.o $(VEGAS_C)
	$(AR) $(ARFLAGS) $(LIB) Vegas.o
	$(RM) Vegas.o

$(LIB)(vegas-f.o): $(VEGAS_F1) $(VEGAS_H)
	$(CC) $(CFLAGS) -c -o vegas-f.o $(VEGAS_F1)
	$(AR) $(ARFLAGS) $(LIB) vegas-f.o
	$(RM) vegas-f.o

$(LIB)(vegas_f.o): $(VEGAS_F2) $(VEGAS_H)
	$(CC) $(CFLAGS) -c -o vegas_f.o $(VEGAS_F2)
	$(AR) $(ARFLAGS) $(LIB) vegas_f.o
	$(RM) vegas_f.o

Vegas: $(VEGAS_MMA) $(VEGAS_SRCS)
	CC=$(CC) mcc $(CFLAGS) -o Vegas $(VEGAS_MMA)
	strip Vegas


SUAVE_C = suave/Suave.c
SUAVE_F1 = suave/suave-f.c
SUAVE_F2 = suave/suave_f.c
SUAVE_MMA = suave/Suave.tm
SUAVE_H = suave/decl.h suave/stddecl.h
SUAVE_SRCS = $(SUAVE_H) suave/util.c suave/debug.c suave/common.c \
  suave/Random.c suave/ChiSquare.c suave/Grid.c suave/Fluct.c \
  suave/Sample.c suave/Integrate.c
 
$(LIB)(Suave.o): $(SUAVE_C) $(SUAVE_SRCS)
	$(CC) $(CFLAGS) -c -o Suave.o $(SUAVE_C)
	$(AR) $(ARFLAGS) $(LIB) Suave.o
	$(RM) Suave.o

$(LIB)(suave-f.o): $(SUAVE_F1) $(SUAVE_H)
	$(CC) $(CFLAGS) -c -o suave-f.o $(SUAVE_F1)
	$(AR) $(ARFLAGS) $(LIB) suave-f.o
	$(RM) suave-f.o

$(LIB)(suave_f.o): $(SUAVE_F2) $(SUAVE_H)
	$(CC) $(CFLAGS) -c -o suave_f.o $(SUAVE_F2)
	$(AR) $(ARFLAGS) $(LIB) suave_f.o
	$(RM) suave_f.o

Suave: $(SUAVE_MMA) $(SUAVE_SRCS)
	CC=$(CC) mcc $(CFLAGS) -o Suave $(SUAVE_MMA)
	strip Suave


DIVONNE_C = divonne/Divonne.c
DIVONNE_F1 = divonne/divonne-f.c
DIVONNE_F2 = divonne/divonne_f.c
DIVONNE_MMA = divonne/Divonne.tm
DIVONNE_H = divonne/decl.h divonne/stddecl.h
DIVONNE_SRCS = $(DIVONNE_H) divonne/util.c divonne/debug.c divonne/common.c \
  divonne/KorobovCoeff.c divonne/Random.c divonne/ChiSquare.c \
  divonne/Rule.c divonne/Sample.c divonne/FindMinimum.c \
  divonne/Explore.c divonne/Split.c divonne/Integrate.c \

$(LIB)(Divonne.o): $(DIVONNE_C) $(DIVONNE_SRCS)
	$(CC) $(CFLAGS) -c -o Divonne.o $(DIVONNE_C)
	$(AR) $(ARFLAGS) $(LIB) Divonne.o
	$(RM) Divonne.o

$(LIB)(divonne-f.o): $(DIVONNE_F1) $(DIVONNE_H)
	$(CC) $(CFLAGS) -c -o divonne-f.o $(DIVONNE_F1)
	$(AR) $(ARFLAGS) $(LIB) divonne-f.o
	$(RM) divonne-f.o

$(LIB)(divonne_f.o): $(DIVONNE_F2) $(DIVONNE_H)
	$(CC) $(CFLAGS) -c -o divonne_f.o $(DIVONNE_F2)
	$(AR) $(ARFLAGS) $(LIB) divonne_f.o
	$(RM) divonne_f.o

Divonne: $(DIVONNE_MMA) $(DIVONNE_SRCS)
	CC=$(CC) mcc $(CFLAGS) -o Divonne $(DIVONNE_MMA)
	strip Divonne


CUHRE_C = cuhre/Cuhre.c
CUHRE_F1 = cuhre/cuhre-f.c
CUHRE_F2 = cuhre/cuhre_f.c
CUHRE_MMA = cuhre/Cuhre.tm
CUHRE_H = cuhre/decl.h cuhre/stddecl.h
CUHRE_SRCS = $(CUHRE_H) cuhre/util.c cuhre/debug.c cuhre/common.c \
  cuhre/ChiSquare.c cuhre/Rule.c cuhre/Integrate.c

$(LIB)(Cuhre.o): $(CUHRE_C) $(CUHRE_SRCS)
	$(CC) $(CFLAGS) -c -o Cuhre.o $(CUHRE_C)
	$(AR) $(ARFLAGS) $(LIB) Cuhre.o
	$(RM) Cuhre.o

$(LIB)(cuhre-f.o): $(CUHRE_F1) $(CUHRE_H)
	$(CC) $(CFLAGS) -c -o cuhre-f.o $(CUHRE_F1)
	$(AR) $(ARFLAGS) $(LIB) cuhre-f.o
	$(RM) cuhre-f.o

$(LIB)(cuhre_f.o): $(CUHRE_F2) $(CUHRE_H)
	$(CC) $(CFLAGS) -c -o cuhre_f.o $(CUHRE_F2)
	$(AR) $(ARFLAGS) $(LIB) cuhre_f.o
	$(RM) cuhre_f.o

Cuhre: $(CUHRE_MMA) $(CUHRE_SRCS)
	CC=$(CC) mcc $(CFLAGS) -o Cuhre $(CUHRE_MMA)
	strip Cuhre


LIBOBJS = \
  Vegas.o vegas-f.o vegas_f.o \
  Suave.o suave-f.o suave_f.o \
  Divonne.o divonne-f.o divonne_f.o \
  Cuhre.o cuhre-f.o cuhre_f.o

$(LIB): $(LIB)($(LIBOBJS))
	-ranlib $(LIB)


demo-fortran: demo-fortran.F $(LIB)
	$(FC) $(FFLAGS) -o demo-fortran demo-fortran.F $(LIB) -lm
	-$(RM) demo-fortran.o

demo-c: demo-c.c cuba.h $(LIB)
	$(CC) $(CFLAGS) -o demo-c demo-c.c $(LIB) -lm


TARFILE = Cuba11.tar.gz

TARDIR = Cuba-1.1

TARCONTENTS = ChangeLog cuba.pdf makefile cuba.h \
  demo-fortran.F demo-c.c demo-c.out.ok demo-math.m testsuite.m \
  $(VEGAS_C) $(VEGAS_F1) $(VEGAS_F2) $(VEGAS_MMA) $(VEGAS_SRCS) \
  $(SUAVE_C) $(SUAVE_F1) $(SUAVE_F2) $(SUAVE_MMA) $(SUAVE_SRCS) \
  $(DIVONNE_C) $(DIVONNE_F1) $(DIVONNE_F2) $(DIVONNE_MMA) $(DIVONNE_SRCS) \
  $(CUHRE_C) $(CUHRE_F1) $(CUHRE_F2) $(CUHRE_MMA) $(CUHRE_SRCS)

tar:
	ln -s . $(TARDIR)
	tar cvfzh $(TARFILE) $(patsubst %,$(TARDIR)/%, $(TARCONTENTS))
	$(RM) $(TARDIR)

pub: tar
	mv -f $(TARFILE) web/
	./mkwebpage

clean:
	-$(RM) $(LIB) $(DEMO) $(MMA) $(LIBOBJS) $(TARFILE) $(TARDIR) demo-c.out

