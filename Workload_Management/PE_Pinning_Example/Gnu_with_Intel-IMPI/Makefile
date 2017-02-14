# Makefile

# Default is set to be for the Intel compilers.  A check is made for the Gnu
# compilers, since the flag to compile OpenMP code differs from Intel's flag.

MPICC=mpicc
CFLAGS=-qopenmp

# Check settings for GNU compilers...

ifeq ($(MPICC_CC),gcc)
  CFLAGS=-fopenmp -openmp-lib compat
endif

.SUFFIXES: .c .f90

all:hybrid_xthi mpi_xthi

hybrid_xthi: xthi.c
	$(MPICC) $(CFLAGS) -o $@ $<

mpi_xthi: mpi_xthi.c
	$(MPICC) -o $@ $<

clean:
	rm -f $(FILES)

