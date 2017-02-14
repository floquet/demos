
all: arrays_write arrays_read

arrays_write: arrays_write.o
	export MPIF90_F90=ifort ; \
   mpif90 -o arrays_write arrays_write.o ${ADIOS_LDFLAGS} ${ADIOS_LIBS}

arrays_write.o: arrays_write.F90
	export MPIF90_F90=ifort ; \
   mpif90 -c ${ADIOS_CPPFLAGS} ${ADIOS_FCFLAGS} arrays_write.F90

arrays_read: arrays_read.o
	export MPIF90_F90=ifort ; \
   mpif90 -o arrays_read arrays_read.o ${ADIOS_LDFLAGS} ${ADIOS_LIBS}

arrays_read.o: arrays_read.F90
	export MPIF90_F90=ifort ; \
   mpif90 -c ${ADIOS_CPPFLAGS} ${ADIOS_FCFLAGS} arrays_read.F90


clean:
	/bin/rm -f *.bp *.o
	/bin/rm -f arrays_write arrays_read

