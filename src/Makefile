
CFLAGS = -Wall -g  -I ../../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -g -cpp 

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = CalculiX.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../SPOOLES.2.2

LIBS = \
       $(DIR)/spooles.a \
	../../LAPACK/lapack-3.12.1/liblapack.a \
	../../LAPACK/lapack-3.12.1/librefblas.a \
	../../ARPACK/libarpack_INTEL.a \
       -lpthread -lm -lc

CalculiX: $(OCCXMAIN) CalculiX.a $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c CalculiX.c; $(FC) -Wall -g -o $@ $(OCCXMAIN) CalculiX.a  $(LIBS) -fopenmp

CalculiX.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
