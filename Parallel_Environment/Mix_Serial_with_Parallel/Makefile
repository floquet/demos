
all: mpi_test.x serial_test.x

mpi_test.x: mpi_test.f90
	mpif90 -o $@ $<

serial_test.x: serial_test.f90
	ifort -o $@ $<


