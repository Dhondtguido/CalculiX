
CFLAGS = -Wall -g  -I /deps/SPOOLES.2.2 -I /deps/SPOOLES.2.2/Misc -I /deps/SPOOLES.2.2/Utilities -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
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

DIR=/deps/SPOOLES.2.2

LIBS = \
       $(DIR)/spooles.a \
	-larpack -lopenblas -llapack -lblas \
       -lpthread -lm -lc

CalculiX: $(OCCXMAIN) CalculiX.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c CalculiX.c; $(FC)  -Wall -g -o $@ $(OCCXMAIN) CalculiX.a $(LIBS) -fopenmp

CalculiX.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

# --- fix: springs.f needs long fixed-form lines ---
# --- fix: springs.f must compile as fixed-form with long lines ---
springs.o: springs.f
	gfortran -Wall -g -cpp -ffixed-form -ffixed-line-length-none -c springs.f -o springs.o
