#!/bin/bash
#
set -v
#
#  For reference, Intel encourages use of their MKL Link Line Advsor at
#  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#
#  The following information is relevant to this example:
#  Topaz runs Parallel Studio XE 2015, Linux OS, no Xeon Phi Co-processors  
#  Intel Fortran compiler, Intel (R) 64 architecture, static linking
#  32-bit integer, sequential execution (ie, no OpenMP),  with SGI MPT
#  communication library

#  The Link Line Advisor yields:

mpif90 -o a.out_scalapack scalapack_exp1.f \
       ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
       -Wl,--start-group \
       ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
       ${MKLROOT}/lib/intel64/libmkl_core.a \
       ${MKLROOT}/lib/intel64/libmkl_sequential.a \
       -Wl,--end-group \
       ${MKLROOT}/lib/intel64/libmkl_blacs_sgimpt_lp64.a \
       -lpthread -lm

# Sometimes an abbreviated form of the compile line works:
# mpif90 -o a.out_scalapack scalapack_exp1.f \
#  -L${MKLROOT}/lib/intel64 \
#  -lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 \
#  -lmkl_core -lmkl_intel_lp64 \
#  -lmkl_lapack95_lp64 -lmkl_sequential

