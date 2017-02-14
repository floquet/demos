


#include <stdio.h>
#include <assert.h>

#include "laio.h"


/*

  -- ScaLAPACK auxiliary routine (version 2.0) --
     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
     and University of California, Berkeley.
     Oct 10, 1996
*/







#if (IDEBUG >= 1)
void lawrite(
#else
static void lawrite(
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

   write a rectangular (m x n) patch  to out of core array A
   starting at A(ia,ja) and copy it from scalapack array B at 
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
    extern void FTN(infog2l) ( fint ia, fint ja, int *descA,
                fint nprow, fint npcol, 
                fint myprow, fint mypcol,
                fint lrindx, fint lcindx,
                fint irow, fint icol );

    extern int lawr(int fd, int irecord, long goffset, int nbytes,
		    void *buf, int currecord,
		    int myid, int nproc);
    extern int lard(int fd, int irecord, long goffset, int nbytes,
		    void *buf, int currecord,
		    int myid, int nproc);

    extern void blacs_pinfo_ (fint myid, fint nproc);
    extern void blacs_gridinfo_ (fint contxt,
		       fint nprow, fint npcol, fint myprow, fint mypcol);
    extern int FTN(numroc) (fint n, fint nb,
			    fint iproc, fint isrcproc, fint nprocs);
    extern int FTN(indxl2g) (fint indxloc,
			fint nb, fint iproc, fint isrcproc, fint nprocs);


    extern int mapblock(int ia, int mmb);
    extern int Findex2(int i, int j, int m, int n);

    extern void *malloc(int);
    extern void free(void *);

    int myid, nproc, nprow, npcol, myprow, mypcol;

    int istart, iend, isize, jstart, jend, jsize;
    int iblock, jblock, iistart, jjstart;
    int irecord, ipos, jpos, ibpos, jbpos;
    int nrow_blocks, ncol_blocks, nbytes, nrow, ncol;
    int row_start, row_end, col_start, col_end;

    int lmmb, lnnb, mb, nb;
    long goffset;

    int mmb, nnb, fd, gm, gn;

    logical isvalid, isinuse, iswrittenout;
    logical hasoverlap,lwriteall, writeall, needrd;

    int ilo,jlo,ihi,jhi, lrindx,lcindx,
	iipos,jjpos, irow,icol,iisize,jjsize, Locp,Locq;

    int nbytes1;

    void *Amat;
    int descA[20] = {-1,-1,-1,-1,-1,
                     -1,-1,-1,-1,-1,
                     -1,-1,-1,-1,-1,
                     -1,-1,-1,-1,-1};

    int Aiodev = (*pAiodev);
    int m = (*pm);
    int n = (*pn);
    int ia = (*pia);
    int ja = (*pja);

    int ib = (*pib);
    int jb = (*pjb);


    blacs_pinfo_ (&myid, &nproc);
    if ( (myid >= 0) && (nproc >= 1)) {
      blacs_gridinfo_ (&(CTXT_A(descA)),
			 &nprow, &npcol, &myprow, &mypcol);
      };


    *info = 0;


    if ((m <= 0) || (n <= 0)) {
	/* nothing to do */
	return;
    };


    isvalid = (1 <= Aiodev) && (Aiodev <= 99);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid Aiodev %d \n", Aiodev);
	};
	*info = -1;
	return;
    };


    isvalid = (1 <= m) && (m <= M_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid m %d M_A(descB) %d ", m, M_A(descB));
	};

	*info = -2;
	return;
    };

    isvalid = (1 <= n) && (n <= N_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid n %d N_A(descB) %d ", n, N_A(descB));
	};

	*info = -3;
	return;
    };

    gm = IOTable[Aiodev].gm;
    gn = IOTable[Aiodev].gn;


    isvalid = (1 <= ia) && (ia + (m - 1) <= gm);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid ia %d %d m gm %d \n", ia, m, gm);
	};

	*info = -4;
	return;
    };


    isvalid = (1 <= ja) && (ja + (n - 1) <= gn);
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid ja %d n %d gn %d \n", ja, n, gn);
	};

	*info = -5;
	return;
    };


    isvalid = (1 <= ib) && (ib + (m - 1) <= M_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid ib %d M_A(descB) %d \n", ib, M_A(descB));
	};

	*info = -7;
	return;
    };


    isvalid = (1 <= jb) && (jb + (n - 1) <= N_A(descB));
    if (!isvalid) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid jb %d N_A(descB) %d \n", jb, N_A(descB));
	};

	*info = -8;
	return;
    };






    /* prepare for main loops */



    isinuse = ((gm >= 1) && (gn >= 1));
    if (!isinuse) {
	if (IDEBUG >= 1) {
	    printf("lawrite: invalid Aiodev %d, not in use \n", Aiodev);
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




    nrow = FTN(numroc) (&(M_A(descA)),
		      &(MB_A(descA)), &myprow, &(RSRC_A(descA)), &nprow);
    ncol = FTN(numroc) (&(N_A(descA)),
		      &(NB_A(descA)), &mypcol, &(CSRC_A(descA)), &npcol);

    if (nrow >= 1) {
	int ii;
	ii = 1;
	row_start = FTN(indxl2g) (&ii, &(MB_A(descA)),
				  &myprow, &(RSRC_A(descA)), &nprow);

	ii = nrow;
	row_end = FTN(indxl2g) (&ii, &(MB_A(descA)),
				&myprow, &(RSRC_A(descA)), &nprow);
    } else {
	row_start = 0;
	row_end = 0;
    };

    if (ncol >= 1) {
	int jj;
	jj = 1;
	col_start = FTN(indxl2g) (&jj, &(NB_A(descA)),
				  &mypcol, &(CSRC_A(descA)), &npcol);

	jj = ncol;
	col_start = FTN(indxl2g) (&jj, &(NB_A(descA)),
				  &mypcol, &(CSRC_A(descA)), &npcol);

    } else {
	col_start = 0;
	col_end = 0;
    };



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
	    writeall = (isize == mmb) && (jsize == nnb);

	    iistart = (iblock - 1) * mmb + 1;
	    jjstart = (jblock - 1) * nnb + 1;

	    ipos = istart - iistart + 1;
	    jpos = jstart - jjstart + 1;

	    ibpos = istart + (ib - ia);
	    jbpos = jstart + (jb - ja);

		/*
		 See if there is any overlap or it entirely
		 covers my local part in
		 A( ipos:ipos+isize-1,  jpos:jpos+jsize-1)
		*/
		FTN( infog2l )( &ipos, &jpos,  descA, &nprow,&npcol,
			&myprow,&mypcol,  &lrindx,&lcindx, &irow,&icol);
	
		ilo = lrindx; jlo = lcindx;
		
		iipos = ipos + isize-1; jjpos = jpos + jsize-1;
                FTN( infog2l )( &iipos, &jjpos,  descA, &nprow,&npcol,
                        &myprow,&mypcol,  &lrindx,&lcindx, &irow,&icol);

		if (irow != myprow) { lrindx = lrindx - 1; };
		if (icol != mypcol) { lcindx = lcindx - 1; };

		ihi = lrindx; jhi = lcindx;
		iisize = MAX(0, (ihi-ilo+1) ); 
		jjsize = MAX(0, (jhi-jlo+1) );

		hasoverlap = ((iisize >= 1) && (jjsize >= 1));
		Locp = FTN(numroc)( &(M_A(descA)), &(MB_A(descA)), 
				&myprow, &(RSRC_A(descA)), &nprow ); 
		Locq = FTN(numroc)( &(N_A(descA)), &(NB_A(descA)),
				&mypcol, &(CSRC_A(descA)), &npcol );
		lwriteall = (Locp == iisize) && (Locq == jjsize);
		iswrittenout = testbit(irecord, IOTable[Aiodev].ondisk);
	
	    if (!writeall) { 
		needrd = (!lwriteall) && hasoverlap && iswrittenout;
		}
	     else { 
		needrd = FALSE; 
		};

/*
	debugging only
*/
	needrd = (!writeall) && iswrittenout;





		if (needrd) {
		    int currecord, newrecord;

		    currecord = IOTable[Aiodev].currecord;
		    newrecord = lard(fd, irecord, goffset, nbytes,
				     Amat, currecord,
				     myid, nproc);
		    IOTable[Aiodev].currecord = newrecord;
		} else if (!lwriteall) {

		    /* fill with zeros */

		    extern void bzero(char *b, int length);

		    bzero((char *) Amat, nbytes);


		};




	    /* perform copy */




	    if (IDEBUG >= 2) {

		if (myid == 0) {
		    printf("lawrite call pdgemr2do \n");
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
		extern void FTN(fpdgemr2d) (fint m, fint n,
				 double *A, fint ia, fint ja, int *descA,
				double *B, fint ib, fint jb, int *descB);

		FTN(fpdgemr2d) (
				   &isize, &jsize,
				   (double *) B, &ibpos, &jbpos, descB,
				   (double *) Amat, &ipos, &jpos, descA);
	    } else if (ctype == 'S') {
		extern void FTN(fpsgemr2d) (fint m, fint n,
				  float *A, fint ia, fint ja, int *descA,
				 float *B, fint ib, fint jb, int *descB);

		FTN(fpsgemr2d) (
				   &isize, &jsize,
				   (float *) B, &ibpos, &jbpos, descB,
				   (float *) Amat, &ipos, &jpos, descA);
	    } else if (ctype == 'C') {
		extern void FTN(fpcgemr2d) (fint m, fint n,
			       complex * A, fint ia, fint ja, int *descA,
			      complex * B, fint ib, fint jb, int *descB);

		FTN(fpcgemr2d) (
				   &isize, &jsize,
				   (complex *) B, &ibpos, &jbpos, descB,
				   (complex *) Amat, &ipos, &jpos, descA);
	    } else if (ctype == 'I') {
		extern void FTN(fpigemr2d) (fint m, fint n,
				    int *A, fint ia, fint ja, int *descA,
				   int *B, fint ib, fint jb, int *descB);

		FTN(fpigemr2d) (
				   &isize, &jsize,
				   (int *) B, &ibpos, &jbpos, descB,
				   (int *) Amat, &ipos, &jpos, descA);
	    } else if (ctype == 'Z') {
		extern void FTN(fpzgemr2d) (fint m, fint n,
			      zcomplex * A, fint ia, fint ja, int *descA,
			     zcomplex * B, fint ib, fint jb, int *descB);

		FTN(fpzgemr2d) (
				   &isize, &jsize,
				   (zcomplex *) B, &ibpos, &jbpos, descB,
				 (zcomplex *) Amat, &ipos, &jpos, descA);
	    } else {
		assert(0);	/* impossible */
	    };


	    /* write back out to disk */

	    {
		int currecord, newrecord;

		currecord = IOTable[Aiodev].currecord;
		newrecord = lawr(fd, irecord, goffset, nbytes,
				 Amat, currecord,
				 myid, nproc);
		IOTable[Aiodev].currecord = newrecord;
	    };


	    /* register this record is now on disk */

	    setbit(irecord, IOTable[Aiodev].ondisk);

	};			/* end for istart */

    };				/* end for jstart */



    free(Amat);


    return;
}








void
 FTN(dlawrite) (
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
    lawrite(pAiodev, pm, pn, pia, pja,
	    B, 'D', pib, pjb, descB, info);
}


void
 FTN(slawrite) (
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
    lawrite(pAiodev, pm, pn, pia, pja,
	    B, 'S', pib, pjb, descB, info);
}


void
 FTN(clawrite) (
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
    lawrite(pAiodev, pm, pn, pia, pja,
	    B, 'C', pib, pjb, descB, info);
}


void
 FTN(zlawrite) (
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
    lawrite(pAiodev, pm, pn, pia, pja,
	    B, 'Z', pib, pjb, descB, info);
}

void
 FTN(ilawrite) (
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
    lawrite(pAiodev, pm, pn, pia, pja,
	    B, 'I', pib, pjb, descB, info);
}
