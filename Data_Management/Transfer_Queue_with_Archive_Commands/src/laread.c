


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




#if (IDEBUG >= 1)
void laread(
#else
static void laread(
#endif
	       fint pAiodev,
	       fint pm,
	       fint pn,
	       fint pia,
	       fint pja,

	       void *B, char ctype,
	       fint pib,
	       fint pjb,
	       int *descB,

	       fint info)
/*
   Purpose:
   ========

   read a rectangular (m x n) patch  from out of core array A
   starting at A(ia,ja) and copy it into scalapack array B at 
   B(ib,jb)


   rely on 'REDIST' to handle the message passing

   assume I/O is much slower than memory to memory copy or message passing


   m,n: intent(in)   patch is (mxn)
   ia,ja: intent(in)  patch start at  A(ia,ja) in out-of-core matrix

   B: intent(out)  scalapack array
   ia,jb: intent(in) patch in array B start at B(ib,jb)

   descB: intent(in) descriptor for array B

   info: intent(out)  0 is ok
 */
{
    extern int lawr(int fd, int irecord, int nbytes,
		    void *buf, int currecord);
    extern int lard(int fd, int irecord, long goffset, int nbytes,
		    void *buf, int currecord,
		    int myid, int nproc);

    extern int mapblock(int ia, int mmb);
    extern int Findex2(int i, int j, int m, int n);

    extern blacs_pinfo_ (fint myid, fint nproc);

    extern void *malloc(int);
    extern void free(void *);



    int myid, nproc;



    int istart, iend, isize, jstart, jend, jsize;
    int iblock, jblock, iistart, jjstart;
    int irecord, ipos, jpos, ibpos, jbpos;
    int nrow_blocks, ncol_blocks, nbytes;

    int nbytes1;
    int lmmb, lnnb, mb, nb;
    long goffset;


    int mmb, nnb, fd, gm, gn;

    logical isvalid, isinuse, iswrittenout;

    void *Amat;
    int descA[20];

    int m, n, ia, ja, ib, jb, Aiodev;


    Aiodev = (*pAiodev);
    m = (*pm);
    n = (*pn);
    ia = (*pia);
    ja = (*pja);

    ib = (*pib);
    jb = (*pjb);




    *info = 0;

    blacs_pinfo_ (&myid, &nproc);


    if ((m <= 0) || (n <= 0)) {
	/* nothing to do */
	return;
    };


    isvalid = (1 <= Aiodev) && (Aiodev <= 99);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid Aiodev %d \n", Aiodev);
	};
	*info = -1;
	return;
    };


    isvalid = (1 <= m) && (m <= M_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid m %d M_A(descB) %d \n", m, M_A(descB));
	};

	*info = -2;
	return;
    };

    isvalid = (1 <= n) && (n <= N_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid n %d N_A(descB) %d \n", n, N_A(descB));
	};

	*info = -3;
	return;
    };

    gm = IOTable[Aiodev].gm;
    gn = IOTable[Aiodev].gn;


    isvalid = (1 <= ia) && (ia + (m - 1) <= gm);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid ia %d gm %d \n", ia, gm);
	};

	*info = -4;
	return;
    };


    isvalid = (1 <= ja) && (ja + (n - 1) <= gn);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid ja %d gn %d \n", ja, gn);
	};

	*info = -5;
	return;
    };






    isvalid = (1 <= ib) && (ib + (m - 1) <= M_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid ib %d m %d M_A(descB) %d \n", 
			ib, m,M_A(descB));
	};

	*info = -7;
	return;
    };


    isvalid = (1 <= jb) && (jb + (n - 1) <= N_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid jb %d N_A(descB) %d \n", jb, N_A(descB));
	};

	*info = -8;
	return;
    };




    /* prepare for main loops */

    gm = IOTable[Aiodev].gm;
    gn = IOTable[Aiodev].gn;

    isinuse = (gm >= 1) && (gn >= 1);
    if (!isinuse) {
	if (IDEBUG >= 1) {
	    printf("laread: invalid Aiodev %d not in use \n", Aiodev);
	};
	*info = -1;
	return;
    };
    mb = IOTable[Aiodev].mb;
    nb = IOTable[Aiodev].nb;

    mmb = IOTable[Aiodev].mmb;
    nnb = IOTable[Aiodev].nnb;
    fd = IOTable[Aiodev].fd;

    lmmb = IOTable[Aiodev].lmmb;
    lnnb = IOTable[Aiodev].lnnb;

    goffset = IOTable[Aiodev].goffset;


    nrow_blocks = mapblock(gm, mmb);
    ncol_blocks = mapblock(gn, nnb);

    if (ctype == 'D') {
	nbytes1 = sizeof(double);
    } else if (ctype == 'S') {
	nbytes1 = sizeof(float);
    } else if (ctype == 'C') {
	nbytes1 = sizeof(complex);
    } else if (ctype == 'Z') {
	nbytes1 = sizeof(zcomplex);
    } else if (ctype == 'I') {
	nbytes1 = sizeof(int);
    } else {
	assert(0);		/* impossible */
    };

    nbytes = lmmb * lnnb * nbytes1;
    Amat = (double *) malloc(nbytes);
    assert(Amat != NULL);

    /* create new descriptor */

    DT_A(descA) = IOTable[Aiodev].dt;
    CTXT_A(descA) = IOTable[Aiodev].ctxt;
    M_A(descA) = mmb;
    N_A(descA) = nnb;
    MB_A(descA) = mb;
    NB_A(descA) = nb;
    RSRC_A(descA) = IOTable[Aiodev].rsrc;
    CSRC_A(descA) = IOTable[Aiodev].csrc;
    LLD_A(descA) = lmmb;




    /* main loop */


    for (jstart = ja; jstart <= ja + n - 1; jstart = jend + 1) {
	jblock = mapblock(jstart, nnb);
	jend = MIN((ja + n - 1), jblock * nnb);
	jsize = jend - jstart + 1;


	for (istart = ia; istart <= ia + m - 1; istart = iend + 1) {

	    iblock = mapblock(istart, mmb);
	    iend = MIN(ia + m - 1, iblock * mmb);
	    isize = iend - istart + 1;

	    irecord = (Findex2(iblock, jblock, nrow_blocks, ncol_blocks) - 1);

	    iswrittenout = testbit(irecord, IOTable[Aiodev].ondisk);
	    if (iswrittenout) {
		int currecord, newrecord;

		currecord = IOTable[Aiodev].currecord;
		newrecord = lard(fd, irecord, goffset, nbytes, Amat, currecord,
				 myid, nproc);
		IOTable[Aiodev].currecord = newrecord;
	    } else {
		/* fill Amat with zeros */

		/*extern void bzero(void *b, int length); */

		bzero((void *) Amat, nbytes);

	    };


	    /* perform copy */

	    iistart = (iblock - 1) * mmb + 1;
	    jjstart = (jblock - 1) * nnb + 1;

	    ipos = istart - iistart + 1;
	    jpos = jstart - jjstart + 1;

	    ibpos = istart + (ib - ia);
	    jbpos = jstart + (jb - ja);



	    if (IDEBUG >= 2) {

		blacs_pinfo_ (&myid, &nproc);

		if (myid == 0) {
		    printf("laread call pdgemr2do \n");
		    printf("istart %d iend %d jstart %d jend %d \n",
			   istart, iend, jstart, jend);
		    printf("iistart %d jjstart %d \n",
			   iistart, jjstart);
		    printf("irecord %d iblock %d jblock %d \n",
			   irecord, iblock, jblock);
		    printf("isize %d jsize %d \n", isize, jsize);
		    printf("ipos %d jpos %d \n", ipos, jpos);
		    printf("ibpos %d jbpos %d \n\n", ibpos, jbpos);
		};
	    };


	    if (ctype == 'D') {
		extern void FTN(fpdgemr2d) (fint isize, fint jsize,
			  double *Amat, fint ipos, fint jpos, int *descA,
			  double *B, fint ibpos, fint jbpos, int *descB);

		FTN(fpdgemr2d) (&isize, &jsize,
				(double *) Amat, &ipos, &jpos, descA,
				(double *) B, &ibpos, &jbpos, descB);

	    } else if (ctype == 'S') {
		extern void FTN(fpsgemr2d) (fint isize, fint jsize,
			   float *Amat, fint ipos, fint jpos, int *descA,
			   float *B, fint ibpos, fint jbpos, int *descB);

		FTN(fpsgemr2d) (&isize, &jsize,
				(float *) Amat, &ipos, &jpos, descA,
				(float *) B, &ibpos, &jbpos, descB);
	    } else if (ctype == 'I') {
		extern void FTN(fpigemr2d) (fint isize, fint jsize,
			     int *Amat, fint ipos, fint jpos, int *descA,
			     int *B, fint ibpos, fint jbpos, int *descB);

		FTN(fpigemr2d) (&isize, &jsize,
				(int *) Amat, &ipos, &jpos, descA,
				(int *) B, &ibpos, &jbpos, descB);
	    } else if (ctype == 'C') {
		extern void FTN(fpcgemr2d) (fint isize, fint jsize,
			complex * Amat, fint ipos, fint jpos, int *descA,
			complex * B, fint ibpos, fint jbpos, int *descB);

		FTN(fpcgemr2d) (&isize, &jsize,
				(complex *) Amat, &ipos, &jpos, descA,
				(complex *) B, &ibpos, &jbpos, descB);
	    } else if (ctype == 'Z') {
		extern void FTN(fpzgemr2d) (fint isize, fint jsize,
		       zcomplex * Amat, fint ipos, fint jpos, int *descA,
		       zcomplex * B, fint ibpos, fint jbpos, int *descB);

		FTN(fpzgemr2d) (&isize, &jsize,
				(zcomplex *) Amat, &ipos, &jpos, descA,
				(zcomplex *) B, &ibpos, &jbpos, descB);
	    } else {
		assert(0);	/* impossible */
	    };






	};			/* end for istart */

    };				/* end for jstart */



    free(Amat);


    return;
}





void
 FTN(dlaread) (
		  fint pAiodev,
		  fint pm,
		  fint pn,
		  fint pia,
		  fint pja,

		  double *B,
		  fint pib,
		  fint pjb,
		  int *descB,

		  fint info) {
    laread(pAiodev, pm, pn, pia, pja,
	   B, 'D', pib, pjb, descB, info);
}


void
 FTN(slaread) (
		  fint pAiodev,
		  fint pm,
		  fint pn,
		  fint pia,
		  fint pja,

		  float *B,
		  fint pib,
		  fint pjb,
		  int *descB,

		  fint info) {
    laread(pAiodev, pm, pn, pia, pja,
	   B, 'S', pib, pjb, descB, info);
}


void
 FTN(claread) (
		  fint pAiodev,
		  fint pm,
		  fint pn,
		  fint pia,
		  fint pja,

		  complex * B,
		  fint pib,
		  fint pjb,
		  int *descB,

		  fint info) {
    laread(pAiodev, pm, pn, pia, pja,
	   B, 'C', pib, pjb, descB, info);
}


void
 FTN(zlaread) (
		  fint pAiodev,
		  fint pm,
		  fint pn,
		  fint pia,
		  fint pja,

		  zcomplex * B,
		  fint pib,
		  fint pjb,
		  int *descB,

		  fint info) {
    laread(pAiodev, pm, pn, pia, pja,
	   B, 'Z', pib, pjb, descB, info);
}

void
 FTN(ilaread) (
		  fint pAiodev,
		  fint pm,
		  fint pn,
		  fint pia,
		  fint pja,

		  double *B,
		  fint pib,
		  fint pjb,
		  int *descB,

		  fint info) {
    laread(pAiodev, pm, pn, pia, pja,
	   B, 'I', pib, pjb, descB, info);
}
