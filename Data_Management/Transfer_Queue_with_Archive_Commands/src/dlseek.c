#include <stdio.h>
#include <assert.h>


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/

/*
 Purpose:
 ========

 32bit integer offset retrict seek access to at most 2gigabytes
 (2^31-1). 'dlseek()' attempts to use extended integers or
 double precision as extended precision (48bits) integers
 to seek past 2 gigabytes.

 This routine is likely to be system dependent.

*/
#ifdef NX

#include <nx.h>
#include <unistd.h>

static void dtoesize(double offset, esize_t * peoffset)
{

    double giga = ((double) 1024) * ((double) 1024) * ((double) 1024);
    double drem;
    long irem, ngiga;

    esize_t eoffset;
    esize_t egiga;
    esize_t eirem;

    ngiga = (long) (offset / giga);
    drem = (offset - ngiga * giga);
    irem = (long) drem;


    egiga = emul(stoe("1073741824"), (long) ngiga);	/* 1024*1024*1024 */
    eirem = emul(stoe("1"), (long) irem);
    eoffset = eadd(egiga, eirem);

    (*peoffset) = eoffset;
}


static void esizetod(esize_t evalue, double *pdvalue)
{
    extern long atol(const char *string);


    double dvalue;
    long lgiga = 1024 * 1024 * 1024;
    long irem, ngiga;

    ngiga = (long) ediv(evalue, (long) lgiga);
    irem = (long) emod(evalue, (long) lgiga);

    dvalue = (double) ngiga;
    dvalue *= (double) lgiga;
    dvalue += (double) irem;

    *pdvalue = dvalue;

}

#endif

double dlseek(int fd, double offset, int whence)
{

/* use double (48bit integer) to get around 2gigabyte limit */

#ifdef NX

    {

	extern esize_t eseek(int fd, esize_t offset, int whence);

	esize_t eoffset, enewoffset;
	double newoffset;

	dtoesize(offset, &eoffset);
	enewoffset = eseek(fd, eoffset, whence);

	esizetod(enewoffset, &newoffset);

	return (newoffset);
    }

#else
    {

	extern long lseek(int fd, long offset, int whence);

	return ((double) lseek(fd, (long) offset, whence));
    }
#endif

}
