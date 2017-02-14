# Settings for sequential binaries and threaded binaries using OpenMP
# Default settings are for the Intel SGI-MPT Programming Environment

CC=mpicc
FF=mpif90

SEROPT=-qno-openmp
OMPOPT=-DOMP -qopenmp

# Intel P.E.
ifeq ($(MPICC_CC),icc)
  SEROPT=-qno-openmp
  OMPOPT=-DOMP -qopenmp
endif

# Gnu P.E.
ifeq ($(MPICC_CC),gcc)
  SEROPT=
  OMPOPT=-DOMP -fopenmp
endif


all: clean corec.x coref.x

openmp: 
	make all "SEROPT=$(OMPOPT)"

corec.x: mainc.o core_info.o
	$(CC) $(SEROPT) -o corec.x mainc.o core_info.o

coref.x: mainf.o core_info.o
	$(FF) $(SEROPT) -o coref.x mainf.o core_info.o

core_info.o: core_info.c
	$(CC) $(SEROPT) -c core_info.c

mainc.o: mainc.c
	$(CC) $(SEROPT) -c mainc.c

mainf.o: mainf.f
	$(FF) $(SEROPT) -c mainf.f

clean:
	rm -f mainc.o mainf.o core_info.o corec.x coref.x
