VERSION = 2.21
.DEFAULT_GOAL := ccx

CFLAGS = -O3 -I/usr/include/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -O3

CC=cc
FC=gfortran

include Makefile.inc

SCCXMAIN = ccx_2.21.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = \
       -lspooles \
	   $(shell pkg-config --libs --static arpack) \
       -lpthread -lm -lc


ccx: $(OCCXMAIN) ccx_$(VERSION).a $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_$(VERSION).c; $(FC) $(FFLAGS) -o ccx_$(VERSION) $(OCCXMAIN) ccx_$(VERSION).a $(LIBS) -fopenmp
	ln -s $(PWD)/ccx_$(VERSION) $(PWD)/ccx

dbg: $(OCCXMAIN) ccx_$(VERSION).a $(LIBS)
	./date.pl; $(CC) -g $(CFLAGS) -c ccx_$(VERSION).c; $(FC) -g $(FFLAGS) -o ccx_$(VERSION) $(OCCXMAIN) ccx_$(VERSION).a $(LIBS) -fopenmp
	ln -s $(PWD)/ccx_$(VERSION) $(PWD)/ccx

ccx_$(VERSION).a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

clean:
	rm -f *.o
	rm -f ccx_$(VERSION) ccx ccx_$(VERSION).a ccx.a
	# make proccess edit this files, with build date
	# but if you change it manually commit the change
	# before any make clean
	git restore ccx_$(VERSION).c ccx_$(VERSION)step.c frd.c
