VERSION = 2.21
.DEFAULT_GOAL := ccx

CFLAGS = -I/${EBROOTSPOOLES}/include/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = 
CC?=cc
FC?=gfortran

include Makefile.inc

SCCXMAIN = ccx_2.21.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = -lspooles -lpthread -larpack -lflexiblas -lm -lc

ccx: $(OCCXMAIN) ccx_$(VERSION).a
        ./date.pl; $(CC) $(CFLAGS) -c ccx_$(VERSION).c; $(FC) $(FFLAGS) -o ccx_$(VERSION) $(OCCXMAIN) ccx_$(VERSION).a $(LIBS) -fopenmp $(LDFLAGS)
        ln -s $(PWD)/ccx_$(VERSION) $(PWD)/ccx

ccx_$(VERSION).a: $(OCCXF) $(OCCXC)
        ar vr $@ $?