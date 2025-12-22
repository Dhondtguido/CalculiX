# ============================
# Bulletproof deps + overrides
# ============================

# Toolchain (overrideable)
CC       ?= cc
FC       ?= gfortran

# Build flags (overrideable)
CPPFLAGS ?=
CFLAGS   ?= -Wall -g
FFLAGS   ?= -Wall -g -cpp
LDFLAGS  ?=
LDLIBS   ?=

# ---- SPOOLES autodetect ----
SPOOLES_DIR_CANDIDATES := /deps/SPOOLES.2.2 ../../SPOOLES.2.2
SPOOLES_DIR ?= $(firstword $(foreach d,$(SPOOLES_DIR_CANDIDATES),$(if $(wildcard $(d)/spooles.a),$(d),)))

ifeq ($(strip $(SPOOLES_DIR)),)
  $(error spooles.a not found. Set SPOOLES_DIR=/path/to/SPOOLES.2.2 (must contain spooles.a))
endif

SPOOLES_LIB   := $(SPOOLES_DIR)/spooles.a
SPOOLES_INCS  := -I $(SPOOLES_DIR) -I $(SPOOLES_DIR)/Misc -I $(SPOOLES_DIR)/Utilities

# ---- ARPACK / BLAS / LAPACK (system libs by default) ----
ARPACK_LIBS ?= -larpack

# Prefer pkg-config openblas if available, fallback otherwise
OPENBLAS_LIBS := $(shell pkg-config --libs openblas 2>/dev/null)
ifeq ($(strip $(OPENBLAS_LIBS)),)
  BLAS_LAPACK_LIBS ?= -lopenblas -llapack -lblas
else
  BLAS_LAPACK_LIBS ?= $(OPENBLAS_LIBS) -llapack -lblas
endif

# Feature macros (keep your existing defines if you have them)
CPPFLAGS += -DARCH=\"Linux\" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT

# Add includes + libs
CPPFLAGS += $(SPOOLES_INCS)
LDLIBS   += $(SPOOLES_LIB) $(ARPACK_LIBS) $(BLAS_LAPACK_LIBS) -lm

# Debug helper
print-vars:
	@echo "CC=$(CC)"
	@echo "FC=$(FC)"
	@echo "SPOOLES_DIR=$(SPOOLES_DIR)"
	@echo "SPOOLES_LIB=$(SPOOLES_LIB)"
	@echo "ARPACK_LIBS=$(ARPACK_LIBS)"
	@echo "BLAS_LAPACK_LIBS=$(BLAS_LAPACK_LIBS)"
	@echo "CPPFLAGS=$(CPPFLAGS)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "LDLIBS=$(LDLIBS)"

# Clean
clean:
	rm -f *.o *.a CalculiX ccx 2>/dev/null || true

CFLAGS += -Wall -g  -I /deps/SPOOLES.2.2 -I /deps/SPOOLES.2.2/Misc -I /deps/SPOOLES.2.2/Utilities -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS += -Wall -g -cpp 

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

springs.o: springs.f
	$(FC) $(FFLAGS) -ffixed-line-length-none -c springs.f -o springs.o

check-deps:
	@test -f "$(SPOOLES_LIB)" || (echo "Missing: $(SPOOLES_LIB)"; exit 1)
	@echo "OK: SPOOLES_LIB=$(SPOOLES_LIB)"
	@$(CC) --version >/dev/null
	@$(FC) --version >/dev/null
