# Makefile

#  Parallel compiler definition and OpenMP flag for Intel compilers

MPICC=mpicc
CFLAGS=-qopenmp

.SUFFIXES: .c .f90

all:hybrid_xthi mpi_xthi

hybrid_xthi: xthi.c
	$(MPICC) $(CFLAGS) -o $@ $<

mpi_xthi: mpi_xthi.c
	$(MPICC) -o $@ $<

clean:
	rm -f $(FILES)

