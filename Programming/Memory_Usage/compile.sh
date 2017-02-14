#! /bin/sh

set -v

export MPICC_CC=icc
MPICC=mpicc
MPIFC=mpif90

/bin/cp measure_memory.c measure_memory_f.c
$MPICC -c -I. -DUNDERSCORECOUNT=1 measure_memory_f.c

/bin/cp measure_memory.c measure_memory_c.c
$MPICC -c -I. -DUNDERSCORECOUNT=0 measure_memory_c.c

$MPICC -c -I. memtest_c.c

$MPIFC -c memtest_f.f90

$MPICC -o memtest_c.exe memtest_c.o measure_memory_c.o

$MPIFC -o memtest_f.exe memtest_f.o measure_memory_f.o

