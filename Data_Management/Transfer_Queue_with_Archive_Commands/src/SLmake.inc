############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   March 20, 1995
#
#  Modified:
#       Jun 4, 2003:    Greatly simplified.  Ripped out all of the "helper"
#                       macros.
#       May 9, 2003:    Reinstated all "factory original" compile flags
#                         (-C and -g).  Doesn't seem to affect output.
#       May 9, 2003:    For BLACSLIB, dropped -lblacssmp, retaining -lblacs
#                         -lblacssmp doesn't seem to do anything for good or ill
#       May 9, 2003:    For BLACSLIB, am adding "-lblacs -lblacssmp" to
#                         the library in BLACSdir.  Kind of klugey, but it
#                         picks up the two last unfound objects: .blacs_pinfo_
#                         and .blacs_grindfo_  (note the trailing underscores
#                         in the names)
#
############################################################################
#
SHELL         = /bin/sh
#
############################################################################
#  PART I.
#  The fortran and C compilers, loaders, and their flags
#
############################################################################
# Supply the names of your preferred FORTRAN and C compilers.
F77           = ifort
CC            = icc

# In provide flags for optimization.
F77FLAGS      = -O3
CCFLAGS       = -O3

# Provide the names of your F77 loader (linker).  If it is the same
# as your compiler, you can leave this macros alone.
F77LOADER     = mpif90 -O3

# Provide any loader flags needed
F77LOADFLAGS  =


#
# DO NOTE EDIT THE USEMPI MACRO!!!!
USEMPI        = -DUsingMpiBlacs


# Choose the correct CDEFS.
# Alpha, SGI, or Sun:
CDEFS         = -DAdd_ $(USEMPI)

# IBM:
#CDEFS         = $(USEMPI)

# Linux:
#CDEFS         = -Df77IsF2C  $(USEMPI)

# T3E:
#CDEFS         = -DUpCase $(USEMPI)


###########################################################################
#  PART II
#  The Libraries
###########################################################################

# Supply the name of your MPI and other normal library, will go at end of loader/linker.

SMPLIB        =
#SMPLIB        = -lmpi



# Supply the name/path of your BLACS FORTRAN init library and you BLACS C init
# library.
BLACSFINIT = 
BLACSCINIT = 



# Supply the name/path of your BLACS library. MPI version and serial

BLACSLIB      =  -lmkl_blacs_sgimpt_lp64

# Supply the name/path of your BLAS library.  (That BLAS, not BLACS!)
# possibilities are provided.  If the blas library is in the LAPACK directory,
# it will probably have a name something like /usr/lib/LAPACK/blas_PWR4.a
BLASLIB       = 


#Supply the name/path of your SCALAPACK library.
SCALAPACKLIB  = -L$(MKLROOT)/lib/intel64/ \
                -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core \
                -lmkl_sequential -lmkl_blacs_sgimpt_lp64 -lpthread -lm


#Supply the name/path of your PBLAS library.  This may be the same as the
# SCALAPACK library above, or it may be a separate library within the SCALAPACK
# directory.
PBLASLIB      =  




#Supply the name/path of your REDIST library.  This may be the same as the SCALAPACK
# library above, or it may be a separate library within the SCALAPACK directory.
REDISTLIB     = 



#Supply the name/path of your SCALAPACK TOOLS library.  This may be the same as the
# SCALAPACK library above, or it may be a separate library within the
# SCALAPACK directory.
TOOLSLIB      = 



#Supply the name/path you want for the optional out-of-core library.  This library
# is produced by "make zcomplex".  It can go anywhere you like.
OUTOFCORELIB = 



# The remaining macros are built from all of the above macros.  If the various
# library names have been supplied correctly, there should be no need to
# edit any of these macros.

CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB)
FBLACSLIB     = $(BLACSFINIT)

PRLIBS        = $(REDISTLIB)  $(CBLACSLIB)
PBLIBS        = $(PBLASLIB) $(TOOLSLIB) $(FBLACSLIB) $(BLASLIB)
LIBS          = $(SCALAPACKLIB)
#LIBS          = $(PBLIBS) $(SCALAPACKLIB) $(PBLIBS)
#LIBS          = $(PBLIBS) $(SCALAPACKLIB)
