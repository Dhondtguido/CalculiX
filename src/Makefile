SPOOLESDIR=../../spooles-2.2

CFLAGS = -Wall -g  -I $(SPOOLESDIR)/include/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -g

CC=gcc
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


LIBS = -L$(NIX_LDFLAGS_FOR_TARGET) \
       $(SPOOLESDIR)/lib/libspooles.a \
	   -larpack -lopenblas \
       -lpthread -lm -lc

CalculiX: $(OCCXMAIN) CalculiX.a 
	$(CC) $(CFLAGS) -c CalculiX.c; $(FC)  -Wall -g -o $@ $(OCCXMAIN) CalculiX.a $(LIBS) -fopenmp

CalculiX.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
