CFLAGS = -Wall -O2 -I/usr/include/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O2

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.21.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = \
       -lspooles \
	   $(shell pkg-config --libs --static arpack) \
       -lpthread -lm -lc

ccx_2.21: $(OCCXMAIN) ccx_2.21.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.21.c; $(FC)  -Wall -O2 -o $@ $(OCCXMAIN) ccx_2.21.a $(LIBS) -fopenmp

ccx_2.21.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
