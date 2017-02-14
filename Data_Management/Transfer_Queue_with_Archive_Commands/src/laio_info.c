#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "laio.h"


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/


void
laio_info_ (
		    Aiodev, mm, nn, mmb, nnb, mb, nb, csrc, rsrc, icontxt)
/*
  Purpose:
  ========

  Return information about out-of-core array.

  These values were set in laopen().

*/

fint Aiodev;
fint mm;
fint nn;
fint mmb;
fint nnb;
fint mb;
fint nb;
fint csrc;
fint rsrc;
fint icontxt;

{

    int iodev;


    iodev = (*Aiodev);
    assert((1 <= iodev) && (iodev <= MAX_IOTABLE));

    *mm = IOTable[iodev].gm;
    *nn = IOTable[iodev].gn;

    *mmb = IOTable[iodev].mmb;
    *nnb = IOTable[iodev].nnb;

    *mb = IOTable[iodev].mb;
    *nb = IOTable[iodev].nb;

    *csrc = IOTable[iodev].csrc;
    *rsrc = IOTable[iodev].rsrc;


    *icontxt = IOTable[iodev].ctxt;

}
