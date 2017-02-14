#include <stdio.h>
#include <assert.h>

#ifndef SEEK_SET
#define SEEK_SET 0
#endif


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/


int lawr(int fd, int irecord, long goffset, 
		int nbytes, void *buf, int currecord,
	 	int myid, int nproc)
/*
 Purpose:
 ========

 Internal auxiliary routine to write a record.

*/

{



    extern int write(int fd, char *buf, int nbyte);
    extern double dlseek(int fd, double offset, int whence);


    int ierr;
    double offset, new_offset;
    int isinterlaced;


    isinterlaced = (goffset < 0);
    if (isinterlaced) {
	/* must perform a seek operation */

	/* recover goffset */
	goffset = -goffset;

	offset = irecord;
	offset *= goffset;
	offset *= nbytes;

	offset += (myid * nbytes);

	new_offset = dlseek(fd, offset, SEEK_SET);
	assert(new_offset == offset);

    } else {

	/* try to avoid an expensive lseek */

	if (irecord != currecord) {
	    offset = (goffset + irecord);
	    offset *= nbytes;
	    new_offset = dlseek(fd, offset, SEEK_SET);

	    assert(new_offset == offset);
	};
    };

    ierr = write(fd, buf, nbytes);
    assert(ierr == nbytes);


    return (irecord + 1);

}
