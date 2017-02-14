
FC = ifort
CC = icc
LIBS = -lmpi
FFLAGS = -g -warn all -check bounds
# FFLAGS = -O3 -warn all

all: mpiio.exe

mpiio.exe: mpiio.f90
	$(FC) $(FFLAGS) -o $@ $< $(LIBS)

clean:
	/bin/rm -f mpiio.exe
	/bin/rm -f *.o
	/bin/rm -f *__genmod.f90  *__genmod.mod
