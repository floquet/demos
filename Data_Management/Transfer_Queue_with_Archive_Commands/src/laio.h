#ifndef LAIO_H
#define LAIO_H 1

#include <stdio.h>
#include <assert.h>

#include "bitvec.h"

#define DT_A( desc ) desc[0]
#define CTXT_A( desc ) desc[1]
#define M_A(desc) desc[2]
#define N_A(desc) desc[3]
#define MB_A(desc) desc[4]
#define NB_A(desc) desc[5]
#define RSRC_A(desc) desc[6]
#define CSRC_A(desc) desc[7]
#define LLD_A(desc) desc[8]

typedef int * fint;
typedef double *fdouble;
typedef float *freal;
typedef int logical;

typedef struct { float real; float imag; } complex;
typedef struct { double real; double imag; } zcomplex;


#ifdef K_AND_R_C

#define name2(a,b) gEnErIc2(a,b)
#define gEnErIc2(a,b) a/**/b

#define name3(a,b,c) gEnErIc3(a,b,c)
#define gEnErIc3(a,b,c) a/**/b/**/c

#define name4(a,b,c,d) gEnErIc4(a,b,c,d)
#define gEnErIc4(a,b,c,d) a/**/b/**/c/**/d

#else /* ANSI_C */

#define name2(a,b) gEnErIc2(a,b)
#define gEnErIc2(a,b) a ## b

#define name3(a,b,c) gEnErIc3(a,b,c)
#define gEnErIc3(a,b,c) a ## b ## c

#define name4(a,b,c,d) gEnErIc4(a,b,c,d)
#define gEnErIc4(a,b,c,d) a ## b ## c ## d

#endif

#ifdef Add_
#define FTN_UNDERSCORES 1
#endif

#ifdef FTN_UNDERSCORES

#define FTN(fname) name2(fname,_)

#else

#define FTN(fname) fname

#endif



#ifndef IDEBUG
#define IDEBUG 1
#endif



/* some data structure */

#define MAX_IOTABLE 100
typedef struct {
	int dt;
	int ctxt;
	int gm;
	int gn;
	int mb;
	int nb;
	int rsrc;
	int csrc;
	int lld;

	int mmb;	/* global template for i/o record */
	int nnb;	/* mmb = lmmb*nprow, nnb=lnnb*npcol */

	int lmmb;	/* local template for i/o record  */
	int lnnb;  
	int fd;

	int currecord;  /* current record */

	char filetype;  /* filetype is 'D' for distributed file,
				       'S' for shared file,
					'I' for interlaced and shared file */

	long goffset;	/* global offset, 0 for distributed file,
			   non-zero for shared files 
			   count in terms of records, not bytes */
	bitvector *ondisk;

	char *filename;
	}  IOTable_t;
extern IOTable_t IOTable[MAX_IOTABLE+1];



#ifndef MOD
#define MOD(ix,n) \
	( (ix) % (n) )
#endif

#ifndef ABS
#define ABS(x) \
	(  ( (x) >= 0 ) ? (x) : (-x) )
#endif

#ifndef MAX
#define MAX(x,y) \
	( ( (x) > (y) ) ? (x) : (y)  )
#endif

#ifndef MIN
#define MIN(x,y) \
    (  ( (x) < (y) ) ? (x) : (y)  )
#endif

#ifndef DIVUP
#define DIVUP(ix,n) \
   ( (MOD(ix,n) == 0) ? ( (ix)/(n) ) : ( (ix)/(n) + 1)  )
#endif



  
#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE (!FALSE)
#endif


#ifdef f77IsF2C
#define laio_info_ laio_info__
#else
#define laio_info_ FTN(laio_info)
#endif


#ifdef f77IsF2C

/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine (which is what the BLACS are written in)
 * for systems where the fortran "compiler" is actually f2c (a fortran
 * to C conversion utility).
 */
/*
 * Initialization routines
 */
#define blacs_pinfo_    blacs_pinfo__
#define blacs_setup_    blacs_setup__
#define blacs_set_      blacs_set__
#define blacs_get_      blacs_get__
#define blacs_gridinit_ blacs_gridinit__
#define blacs_gridmap_  blacs_gridmap__
/*
 * Destruction routines
 */
#define blacs_freebuff_ blacs_freebuff__
#define blacs_gridexit_ blacs_gridexit__
#define blacs_abort_    blacs_abort__
#define blacs_exit_     blacs_exit__
/*
 * Informational & misc.
 */
#define blacs_gridinfo_ blacs_gridinfo__
#define blacs_pnum_     blacs_pnum__
#define blacs_pcoord_   blacs_pcoord__
#define blacs_barrier_  blacs_barrier__

#endif



#endif
