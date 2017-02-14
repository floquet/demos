
all: adios_global

adios_global: adios_global.o
	export MPIF90_F90=ifort ; \
   mpif90 -o adios_global adios_global.o ${ADIOS_LDFLAGS} ${ADIOS_LIBS}

adios_global.o: adios_global.F90
	export MPIF90_F90=ifort ; \
   mpif90 -c ${ADIOS_CPPFLAGS} ${ADIOS_FCFLAGS} adios_global.F90


clean:
	/bin/rm -f *.bp *.o
	/bin/rm -f adios_global
	touch adios_global.bp
