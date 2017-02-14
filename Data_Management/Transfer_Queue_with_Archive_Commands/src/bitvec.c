
#include <stdio.h>
#include <assert.h>

#include "bitvec.h"

#define nbitsperbyte 8
#define nbitspervector (sizeof(bitvector)*nbitsperbyte)

#define mask( jthbit ) \
	((bitvector) ( ( (bitvector) 1)  << (jthbit) ))


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/

/*
 Purpose:
 ========

 Routines to manipulate a bit vector.

*/


bitvector *newbitvec(int nbits)
{

    extern void *calloc(int, int);

    int nword;

    nword = (nbits / nbitspervector) + 2;	/* add two just to be safe */

    /* initialy all bits are zero */
    return ((bitvector *) calloc(1, sizeof(bitvector) * nword));


}

void freebitvec(bitvector * bitvec)
{
    extern void free(void *);

    assert(bitvec != NULL);
    free((void *) bitvec);
}

void setbit(int ithbit, bitvector bitvec[])
{
    int nword, jthbit;

    assert(0 <= ithbit);

    nword = ithbit / nbitspervector;

    jthbit = ithbit - nword * nbitspervector;

    assert((0 <= jthbit) & (jthbit < nbitspervector));

    bitvec[nword] |= mask(jthbit);
}



void clrbit(int ithbit, bitvector bitvec[])
{
    int nword, jthbit;

    assert(0 <= ithbit);

    nword = ithbit / nbitspervector;

    jthbit = ithbit - nword * nbitspervector;

    assert((0 <= jthbit) & (jthbit < nbitspervector));

    bitvec[nword] &= ~(mask(jthbit));
}

int testbit(int ithbit, bitvector bitvec[])
{

    int nword, jthbit;

    assert(0 <= ithbit);

    nword = ithbit / nbitspervector;

    jthbit = ithbit - nword * nbitspervector;

    assert((0 <= jthbit) & (jthbit < nbitspervector));


    return ((bitvec[nword] & mask(jthbit)) >> jthbit);
}
