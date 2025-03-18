PASTIXDIR=../../pastix-6.4.0

#CFLAGS = -Wall -g $(NIX_CFLAGS_COMPILE) -I $(SPOOLESDIR)/include/spooles  -I $(PASTIXDIR)/include -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DPASTIX
CFLAGS = -Wall -O0 -g -fopenmp $(NIX_CFLAGS_COMPILE) -DARCH="Linux" -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DPASTIX -DINTSIZE64
FFLAGS = -Wall -O0 -g -fdefault-integer-8

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


LIBS = $(NIX_LDFLAGS_FOR_TARGET) \
	   -lpastix -lpastix_starpu -lpastix -lpastix_kernels -lspm -lstarpu-1.4 -lhwloc -lscotch -lscotcherrexit \
	   -larpack -lblas \
       -lpthread -lm -lc
	   #-L$(PASTIXDIR)/lib -lpastix -lpastix_kernels -lspm -lhwloc -lscotch -lscotcherrexit \

CalculiX: $(OCCXMAIN) CalculiX.a 
	$(CC) $(CFLAGS) -c CalculiX.c; $(FC)  -Wall -g -o $@ $(OCCXMAIN) CalculiX.a $(LIBS) -fopenmp

CalculiX.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm -f *.o *.a
