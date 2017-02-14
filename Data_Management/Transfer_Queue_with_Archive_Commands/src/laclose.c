#include <stdio.h>
#include <assert.h>
#include "laio.h"


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/


extern void free(void *);
extern int unlink(char *path);
extern int close(int fd);

void
 FTN(laclose) (
		  fint piodev,
		  char *pkeep,
		  fint pmyid,
		  fint pnproc,
		  fint info
)
/*
   Purpose:
   ========

   Close a file associated with an out-of-core matrix.

   iodev: intent(in)
   keep: 'K' to keep the file, otherwise,
	  the file will be unlink'ed after closing

   
 */

{

    char *filename;
    logical isvalid, isinuse;
    int ierr, fd;

    int gm, gn;

    int iodev;
    char ckeep;
    logical iskeep;

    logical isshared, isdistributed, isinterlaced;
    char filetype;

    int myid, nproc;



    *info = 0;

    myid = (*pmyid);
    nproc = (*pnproc);

    iodev = (*piodev);
    isvalid = (1 <= iodev) && (iodev <= 99);
    if (!isvalid) {
	*info = -1;
	return;
    };

    gm = IOTable[iodev].gm;
    gn = IOTable[iodev].gn;

    isinuse = ((gm >= 1) && (gn >= 1));

    if (!isinuse) {
	*info = -1;
	return;
    };


    ckeep = (*pkeep);
    iskeep = (ckeep == 'K') || (ckeep == 'k');



    /* clear out data structure */


    IOTable[iodev].dt = 0;
    IOTable[iodev].ctxt = 0;
    IOTable[iodev].gm = 0;
    IOTable[iodev].gn = 0;
    IOTable[iodev].mb = 0;
    IOTable[iodev].nb = 0;

    IOTable[iodev].mmb = 0;
    IOTable[iodev].nnb = 0;
    IOTable[iodev].lmmb = 0;
    IOTable[iodev].lnnb = 0;

    IOTable[iodev].rsrc = 0;
    IOTable[iodev].csrc = 0;
    IOTable[iodev].lld = 0;


    filename = IOTable[iodev].filename;

    fd = IOTable[iodev].fd;
    ierr = close(fd);
    isvalid = (ierr == 0);
    if (!isvalid) {

	if (IDEBUG >= 1) {
	    printf("laclose: error in closing file %s\n",
		   filename);
	};
	*info = -1;
	return;
    };


    IOTable[iodev].fd = 0;
    IOTable[iodev].currecord = -1;

    filetype = IOTable[iodev].filetype;
    isshared = (filetype == 'S') || (filetype == 's');
    isdistributed = (filetype == 'D') || (filetype == 'd');
    isinterlaced = (filetype == 'I') || (filetype == 'i');

    if (!iskeep) {
	/* don't keep this large file */

	if (isshared || isinterlaced) {
	    /* every processor shares only one file */

	    if (myid == 0) {
		ierr = unlink(filename);
		assert(ierr == 0);
	    };
	} else {
	    /* file is distributed, every processor has a unique filename */

	    ierr = unlink(filename);
	    assert(ierr == 0);
	};
    };


    /* release memory */

    freebitvec(IOTable[iodev].ondisk);
    IOTable[iodev].ondisk = NULL;

    free(filename);
    IOTable[iodev].filename = NULL;

    IOTable[iodev].goffset = 0;
    IOTable[iodev].filetype = ' ';

    return;
}
