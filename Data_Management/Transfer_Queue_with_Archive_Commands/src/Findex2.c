#include <stdio.h>
#include <assert.h>

/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/

int Findex2(int i, int j, int m, int n)
{

/*
 Purpose:
 ========

 Generate offset into a matrix given fortran indexing.

*/
    assert((1 <= i) && (i <= m));
    assert((1 <= j) && (j <= n));

    return (i + (j - 1) * m);
}
