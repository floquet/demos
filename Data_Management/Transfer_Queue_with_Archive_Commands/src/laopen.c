#include <stdio.h>
#include <assert.h>
#include <string.h>

/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/



#include <fcntl.h>


#if defined(__PARAGON__) 
/* Most likely this is an NX system */
#define NX 1
#endif



#include "laio.h"

IOTable_t IOTable[MAX_IOTABLE+1];



/* standard unix open */
#define GOPEN( fname, flags, mode ) \
	open(fname,flags,mode)

#define SETSIZE(fd, nbytes )


#ifdef NX
#define FILE_PREFIX ""
#include "nx.h"

#ifdef M_ASYNC
#define M_TYPE M_ASYNC
#else
#define M_TYPE M_UNIX
#endif

/* special gopen for NX */

#undef GOPEN
#define GOPEN( fname, flags, mode ) \
	gopen( fname, flags, M_TYPE, mode )


#undef SETSIZE
#define SETSIZE( fd, nbytes ) { \
	if (myid == 0) { \
	  (void) lsize( fd, (long) nbytes, (long) 0 ); \
	  }; \
 	}
#endif

#ifdef SP2
#define FILE_PREFIX ""
#endif

#ifndef FILE_PREFIX
#define FILE_PREFIX ""
#endif

/*
   open a file  for storing an out of core  array 
 */

#ifndef S_IREAD
#define S_IREAD         0000400	/* read permission, owner */
#endif

#ifndef S_IWRITE
#define S_IWRITE        0000200	/* write permission, owner */
#endif

#ifndef S_IRUSR
#define         S_IRUSR 0000400	/* read permission, owner */
#endif

#ifndef S_IWUSR
#define         S_IWUSR 0000200	/* write permission, owner */
#endif

#ifndef S_IRGRP
#define         S_IRGRP 0000040	/* read permission, group */
#endif


#ifndef S_IWGRP
#define         S_IWGRP 0000020	/* write permission, grougroup */
#endif

#ifndef S_IROTH
#define         S_IROTH 0000004	/* read permission, other */
#endif

#ifndef S_IWOTH
#define         S_IWOTH 0000002	/* write permission, other */
#endif


#ifndef O_RDWR
#define O_RDWR 2
#endif

#ifndef O_CREAT
#define O_CREAT  0x0200
#endif


void
 FTN(laopen) (
		 fint piodev,
		 char *filetype,
		 char *filename,
		 int *desc,
		 fint plmmb,
		 fint plnnb,
		 fint pmmb,
		 fint pnnb,
		 fint pmyid,
		 fint pnproc,
		 fint info)
#define iodevpos0 1
#define filetypepos0 2
#define filenamepos0 3
#define descpos0 4
#define lmmbpos0 5
#define lnnbpos0 6
#define mmbpos0 7
#define nnbpos0 8
#define myidpos0 9
#define nprocpos0 10
/* 
   Purpose:
   ========

   Open a file for operations with out-of-core matrix.


   iodev: intent(in)     similar to fortran io unit  1..99
   filetype: intent(in)  
		'S' for shared file, 
		'D' for distributed file
		'I' for shared 'interleaved' file.

   filename: intent(in)  'name of file'

		If filename(1:1) begins with a '/',
		then it assumes an absolute path;

		Otherwise, a suffix path is added.
		On the Paragon, it may be '/pfs', on
		the SP2, it may be '/tmp'.


   desc:  intent(in) descriptor for scalapack array

   lmmb,lnnb: intent(in),  local shape of block stored on disk
		lmmb should be a multiple of mb,
		lnnb should be a multiple of nb

   mmb, nnb : intent(in),   global shape of block stored on disk
   mmb = lmmb*nprow, nnb = lnnb*npcol

   myid : intent(in),  my processor id  (0..nproc-1)
   nproc: intent(in),  number of processors involved in context

   info: intent(out)   0 is ok

 */
{

    extern char *strdup(const char *);
    /* extern int open( const char *path, int flags, mode_t mode ); */

#ifndef SEEK_SET
#define SEEK_SET 0
#endif


    extern int mapblock(int ia, int mmb);

#define MAXLENGTH 256


    int iodev = (*piodev);

    int lmmb = (*plmmb);
    int lnnb = (*plnnb);

    int mmb = (*pmmb);
    int nnb = (*pnnb);

    int myid = (*pmyid);
    int nproc = (*pnproc);

    char ctype = *filetype;

    int m, n,  mb, nb;

    char fname[MAXLENGTH + 1];

    logical isvalid, isshared, isdistributed, isinterlaced, spaceok,
     isabsolute_path;
    logical isavailable;

    long nbytes,goffset;
    int fd, flags, mode, nrecords, nrow_blocks, ncol_blocks;

    /* check parameters */

    *info = 0;

    isshared = (ctype == 'S') || (ctype == 's');
    isdistributed = (ctype == 'D') || (ctype == 'd');
    isinterlaced = (ctype == 'I') || (ctype == 'i');

    isvalid = (isshared || isdistributed || isinterlaced);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: invalid ctype %c \n", ctype);
	};

	*info = -filetypepos0;
	return;
    };


    m = M_A(desc);
    n = N_A(desc);
    mb = MB_A(desc);
    nb = NB_A(desc);

    isvalid = ((1 <= (mmb)) && ((mmb) <= (M_A(desc))));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: warning ** mmb %d m %d \n",
		   mmb, M_A(desc));
	};
	*info = -mmbpos0;
    };

    isvalid = (1 <= (nnb)) && ((nnb) <= (N_A(desc)));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: warning ** nnb %d n %d \n",
		   nnb, N_A(desc));
	};
	*info = -nnbpos0;
    };

    isvalid = (1 <= (nproc));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: invalid nproc %d \n", nproc);
	};
	*info = -nprocpos0;
	return;
    };

    isvalid = ((0 <= (myid)) && ((myid) < (nproc)));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: invalid myid %d \n", myid);
	};
	*info = -myidpos0;
	return;
    };

    /* modify file name */

    isabsolute_path = (filename[0] == '/');
    if (isabsolute_path) {
	strncpy(fname, filename, MAXLENGTH);
    } else {
	/* add  a prefix */
	strcpy(fname, FILE_PREFIX);
	strncat(fname, filename, MAXLENGTH - strlen(FILE_PREFIX));
    }



    if (isdistributed) {

	/* each processor has a unique file name */

	spaceok = ((strlen(fname) + strlen(".0123")) <= MAXLENGTH);
	if (!spaceok) {
	    if (IDEBUG >= 1) {
		printf("laopen: filename %s too long \n", fname);
	    };
	    *info = -filenamepos0;
	    return;
	};

	sprintf(&(fname[strlen(fname)]), ".%04d", (myid));
    };



    /* attempt to open file with default mode */


    flags = (O_RDWR | O_CREAT);
    mode = S_IREAD | S_IWRITE | S_IRUSR | S_IWUSR |
	S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;

    if (isshared || isinterlaced) {

	fd = GOPEN(fname, flags, mode);

	nbytes = MAX(0, mmb*nnb*sizeof(double) );
	SETSIZE( fd,  nbytes );
	    

    } else {

	fd = open(fname, flags, mode);

    };


    isvalid = (fd >= 1);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: error in opening file '%s'\n", fname);
	};

	*info = -filenamepos0;
	return;
    };


    nrow_blocks = mapblock(m, mmb);
    ncol_blocks = mapblock(n, nnb);

    if (isshared) {

	/* goffset in terms of records, not bytes */

	goffset = nrow_blocks * ncol_blocks;
	goffset *= myid;
    } else if (isinterlaced) {

	/* store a negative value to indicate this is an interlaced file  */

	goffset = nrow_blocks * ncol_blocks;
	goffset = -goffset;
    } else {
	goffset = 0;
    }



    /* setup internal data structures */



    isvalid = (1 <= iodev) && (iodev <= 99);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: invalid iodev %d \n", iodev);
	};
	*info = -iodevpos0;
	return;
    };


    isavailable = (IOTable[iodev].gm == 0) &&
	(IOTable[iodev].gn == 0);

    if (!isavailable) {
	if (IDEBUG >= 1) {
	    printf("laopen: invalid iodev %d alread in use \n", iodev);
	    printf("IOTable[iodev].gm %d  IOTable[iodev].gn %d \n",
		IOTable[iodev].gm, IOTable[iodev].gn );
	};
	*info = -iodevpos0;
	return;
    };

    /* double check */

    isvalid = (MOD(lmmb, mb) == 0) &&
	(MOD(lnnb, nb) == 0) &&
	(MOD(mmb, lmmb) == 0) &&
	(MOD(nnb, lnnb) == 0);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laopen: mb %d nb %d lmmb %d lnn %d mmb %d nnb %d \n",
		   mb, nb, lmmb, lnnb, mmb, nnb);
	};
	if (MOD(lmmb,mb) != 0) { *info = -lmmbpos0; }
	else if (MOD(lnnb,nb) != 0) { *info = -lnnbpos0; }
	else if (MOD(mmb,lmmb) != 0) { *info = -mmbpos0; }
	else if (MOD(nnb,lnnb) != 0) { *info = -nnbpos0; };

	return;
    };





    IOTable[iodev].dt = DT_A(desc);
    IOTable[iodev].ctxt = CTXT_A(desc);
    IOTable[iodev].gm = M_A(desc);
    IOTable[iodev].gn = N_A(desc);
    IOTable[iodev].mb = MB_A(desc);
    IOTable[iodev].nb = NB_A(desc);
    IOTable[iodev].rsrc = RSRC_A(desc);
    IOTable[iodev].csrc = CSRC_A(desc);
    IOTable[iodev].lld = LLD_A(desc);


    IOTable[iodev].lmmb = lmmb;
    IOTable[iodev].lnnb = lnnb;

    IOTable[iodev].mmb = mmb;
    IOTable[iodev].nnb = nnb;


    IOTable[iodev].filetype = (*filetype);
    IOTable[iodev].goffset = goffset;
    IOTable[iodev].fd = fd;

    IOTable[iodev].currecord = -1;

    nrecords = (nrow_blocks * ncol_blocks + 1);		/* add one just to be safe */
    IOTable[iodev].ondisk = newbitvec(nrecords);
    IOTable[iodev].filename = strdup((&(fname[0])));




    return;
}
