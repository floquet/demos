#include <stdio.h>
#include <assert.h>

#include "laio.h"


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/


int mapblock(int ia, int mmb)
/*
 Purpose:
 ========

 Internal auxiliary routine to compute block offset.

*/
{
    int ival;
    logical isvalid;

    ival = DIVUP(ia, mmb);

    isvalid = ((ival - 1) * mmb < ia) && (ia <= ival * mmb);

    if (!isvalid) {
	printf("mapblock: ia %d mmb %d ival %d \n",
	       ia, ival, mmb);
    };
    assert(isvalid);

    return (ival);
}
