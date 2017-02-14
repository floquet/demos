#
# use FTYPE=for for final distribution
#
FTYPE=for

PMATGEN_FILES=\
	pcbmatgen.$(FTYPE) \
	pcmatgen.$(FTYPE) \
	pdbmatgen.$(FTYPE) \
	pdmatgen.$(FTYPE) \
	pmatgeninc.$(FTYPE) \
	psbmatgen.$(FTYPE) \
	psmatgen.$(FTYPE) \
	pzbmatgen.$(FTYPE)  \
	pzmatgen.$(FTYPE) 




# GEMM_FILES= pzgemm_.c psgemm_.c pdgemm_.c pcgemm_.c
GEMM_FILES= 



SHIFT_FILES=\
	pishift.$(FTYPE) \
	pdshift.$(FTYPE) psshift.$(FTYPE) \
	pzshift.$(FTYPE) pcshift.$(FTYPE)

AUXILIARY_FILES=\
	$(GEMM_FILES) \
	$(SHIFT_FILES) \
	$(PMATGEN_FILES) \
	setupblacs.$(FTYPE) \
	fpcgemr2d.$(FTYPE) fpdgemr2d.$(FTYPE) \
	fpsgemr2d.$(FTYPE) fpzgemr2d.$(FTYPE) fpigemr2d.$(FTYPE) \
	checkalign.$(FTYPE) \
	assert.$(FTYPE) \
	ifill.$(FTYPE) icopy.$(FTYPE) \
	dfill.$(FTYPE) sfill.$(FTYPE) cfill.$(FTYPE) zfill.$(FTYPE) \
	\
	numrocinv.$(FTYPE) numroc2.$(FTYPE) \
	chcalsize.$(FTYPE) \
	calcolsize.$(FTYPE) \
	pfdescinit.$(FTYPE) \
	pfmaxsize.$(FTYPE) pfmax2size.$(FTYPE) 

ZCOMPLEX_FILES=\
	pgztrsm.$(FTYPE) \
	pfzmatgen.$(FTYPE) pfzlaprnt.$(FTYPE) pfzlaprnt2.$(FTYPE) pfzgemm.$(FTYPE) \
	pfztrsm.$(FTYPE) pfzlatrsm.$(FTYPE)  \
	\
	pfzgeqrf.$(FTYPE) pfzqrfact2.$(FTYPE) pfzqrapply.$(FTYPE) \
	pfzgeqrs.$(FTYPE) pfzqrsolve.$(FTYPE) pfzunmqr.$(FTYPE) \
	pfzcopytau.$(FTYPE) \
	\
	pfzgetrf.$(FTYPE) pfzgetf2.$(FTYPE) pfzoclu.$(FTYPE) \
	pfzgelup.$(FTYPE) pfzgetrs.$(FTYPE) pfzlusolve.$(FTYPE) \
	\
	pfzpotrf.$(FTYPE) pfzpotrs.$(FTYPE) \
	pfzchfact.$(FTYPE) \
	pfzuchfact.$(FTYPE) pfzlchfact.$(FTYPE)


COMPLEX_FILES=\
	pgctrsm.$(FTYPE) \
	pfcmatgen.$(FTYPE) pfclaprnt.$(FTYPE) pfclaprnt2.$(FTYPE) pfcgemm.$(FTYPE) \
	pfctrsm.$(FTYPE) pfclatrsm.$(FTYPE)  \
	\
	pfcgeqrf.$(FTYPE) pfcqrfact2.$(FTYPE) pfcqrapply.$(FTYPE) \
	pfcgeqrs.$(FTYPE) pfcqrsolve.$(FTYPE) pfcunmqr.$(FTYPE) \
	pfccopytau.$(FTYPE) \
	\
	pfcgetrf.$(FTYPE) pfcgetf2.$(FTYPE)  pfcoclu.$(FTYPE) \
	pfcgelup.$(FTYPE) pfcgetrs.$(FTYPE) pfclusolve.$(FTYPE) \
	\
	pfcpotrf.$(FTYPE) pfcpotrs.$(FTYPE) \
	pfcchfact.$(FTYPE) \
	pfcuchfact.$(FTYPE) pfclchfact.$(FTYPE)


DOUBLE_FILES=\
	pgdtrsm.$(FTYPE) \
	pfdmatgen.$(FTYPE) pfdlaprnt.$(FTYPE) pfdlaprnt2.$(FTYPE) pfdgemm.$(FTYPE) \
	pfdtrsm.$(FTYPE) pfdlatrsm.$(FTYPE)  \
	\
	pfdgeqrf.$(FTYPE) pfdqrfact2.$(FTYPE) pfdqrapply.$(FTYPE) \
	pfdgeqrs.$(FTYPE) pfdqrsolve.$(FTYPE) pfdormqr.$(FTYPE) \
	pfdcopytau.$(FTYPE) \
	\
	pfdgetrf.$(FTYPE) pfdgetf2.$(FTYPE)  pfdoclu.$(FTYPE) \
	pfdgelup.$(FTYPE) pfdgetrs.$(FTYPE) pfdlusolve.$(FTYPE) \
	\
	pfdpotrf.$(FTYPE) pfdpotrs.$(FTYPE) \
	pfdchfact.$(FTYPE) \
	pfduchfact.$(FTYPE) pfdlchfact.$(FTYPE)

REAL_FILES=\
	pgstrsm.$(FTYPE) \
	pfsmatgen.$(FTYPE) pfslaprnt.$(FTYPE) pfslaprnt2.$(FTYPE) pfsgemm.$(FTYPE) \
	pfstrsm.$(FTYPE) pfslatrsm.$(FTYPE)  \
	\
	pfsgeqrf.$(FTYPE) pfsqrfact2.$(FTYPE) pfsqrapply.$(FTYPE) \
	pfsgeqrs.$(FTYPE) pfsqrsolve.$(FTYPE) pfsormqr.$(FTYPE) \
	pfscopytau.$(FTYPE) \
	\
	pfsgetrf.$(FTYPE) pfsgetf2.$(FTYPE) pfsoclu.$(FTYPE) \
	pfsgelup.$(FTYPE) pfsgetrs.$(FTYPE) pfslusolve.$(FTYPE) \
	\
	pfspotrf.$(FTYPE) pfspotrs.$(FTYPE) \
	pfschfact.$(FTYPE) \
	pfsuchfact.$(FTYPE) pfslchfact.$(FTYPE)

DRIVER_FILES=\
	testddriver.$(FTYPE) testsdriver.$(FTYPE) \
	testzdriver-lrg.$(FTYPE) testcdriver.$(FTYPE) \
	testzdriver-std.$(FTYPE) testcdriver.$(FTYPE) \
	driverinfo.$(FTYPE) testdriver.in


# GEMM_OBJECTS= pzgemm_.o psgemm_.o pdgemm_.o pcgemm_.o
GEMM_OBJECTS= 

PMATGEN_OBJECTS= \
        pcbmatgen.o \
        pcmatgen.o \
        pdbmatgen.o \
        pdmatgen.o \
        pmatgeninc.o \
        psbmatgen.o \
        psmatgen.o \
        pzbmatgen.o \
        pzmatgen.o 

SHIFT_OBJECTS=\
	pishift.o \
	pdshift.o psshift.o \
	pzshift.o pcshift.o


HEADER_FILES=bitvec.h pblas.h laio.h namecat.h
IOAUX_FILES=mapblock.c lard.c lawr.c bitvec.c dlseek.c Findex2.c
IO_FILES=\
	$(HEADER_FILES) \
	$(IOAUX_FILES) \
	laread.c lawrite.c  laopen.c laclose.c laio_info.c 

DISTRIBUTION_FILES=\
	$(AUXILIARY_FILES) \
	$(IO_FILES) \
	$(DRIVER_FILES) \
	$(ZCOMPLEX_FILES) $(COMPLEX_FILES) $(DOUBLE_FILES) $(REAL_FILES)


#
#
AUXILIARY_OBJECTS=\
	$(GEMM_OBJECTS) \
	$(SHIFT_OBJECTS) \
	$(PMATGEN_OBJECTS) \
	setupblacs.o \
        fpcgemr2d.o fpdgemr2d.o \
	fpsgemr2d.o fpzgemr2d.o fpigemr2d.o \
        checkalign.o \
        assert.o \
	ifill.o icopy.o \
	dfill.o sfill.o cfill.o zfill.o \
        \
        numrocinv.o numroc2.o \
	chcalsize.o \
	calcolsize.o \
	pfdescinit.o \
	pfmaxsize.o pfmax2size.o

ZCOMPLEX_OBJECTS=\
	pgztrsm.o \
	pfzmatgen.o pfzlaprnt.o pfzlaprnt2.o pfzgemm.o \
	pfztrsm.o pfzlatrsm.o  \
	\
	pfzgeqrf.o pfzqrfact2.o pfzqrapply.o \
	pfzgeqrs.o pfzqrsolve.o pfzunmqr.o \
	pfzcopytau.o \
	\
	pfzgetrf.o pfzgetf2.o pfzoclu.o \
	pfzgelup.o pfzgetrs.o pfzlusolve.o \
	\
	pfzpotrf.o pfzpotrs.o \
	pfzchfact.o \
	pfzuchfact.o pfzlchfact.o


COMPLEX_OBJECTS=\
	pgctrsm.o \
	pfcmatgen.o pfclaprnt.o pfclaprnt2.o pfcgemm.o \
	pfctrsm.o pfclatrsm.o  \
	\
	pfcgeqrf.o pfcqrfact2.o pfcqrapply.o \
	pfcgeqrs.o pfcqrsolve.o pfcunmqr.o \
	pfccopytau.o \
	\
	pfcgetrf.o pfcgetf2.o pfcoclu.o \
	pfcgelup.o pfcgetrs.o pfclusolve.o \
	\
	pfcpotrf.o pfcpotrs.o \
	pfcchfact.o \
	pfcuchfact.o pfclchfact.o


DOUBLE_OBJECTS=\
	pgdtrsm.o \
	pfdmatgen.o pfdlaprnt.o pfdlaprnt2.o pfdgemm.o \
	pfdtrsm.o pfdlatrsm.o  \
	\
	pfdgeqrf.o pfdqrfact2.o pfdqrapply.o \
	pfdgeqrs.o pfdqrsolve.o pfdormqr.o \
	pfdcopytau.o \
	\
	pfdgetrf.o pfdgetf2.o pfdoclu.o \
	pfdgelup.o pfdgetrs.o pfdlusolve.o \
	\
	pfdpotrf.o pfdpotrs.o \
	pfdchfact.o \
	pfduchfact.o pfdlchfact.o


REAL_OBJECTS=\
	pgstrsm.o \
	pfsmatgen.o pfslaprnt.o pfslaprnt2.o pfsgemm.o \
	pfstrsm.o pfslatrsm.o  \
	\
	pfsgeqrf.o pfsqrfact2.o pfsqrapply.o \
	pfsgeqrs.o pfsqrsolve.o pfsormqr.o \
	pfscopytau.o \
	\
	pfsgetrf.o pfsgetf2.o pfsoclu.o \
	pfsgelup.o pfsgetrs.o pfslusolve.o \
	\
	pfspotrf.o pfspotrs.o \
	pfschfact.o \
	pfsuchfact.o pfslchfact.o
DRIVER_OBJECTS=\
	testzdriver.o \
	testcdriver.o \
	testddriver.o \
	testsdriver.o \
	driverinfo.o 

IOAUX_OBJECTS=mapblock.o lard.o lawr.o bitvec.o dlseek.o Findex2.o
IO_OBJECTS=\
	$(IOAUX_OBJECTS) \
	laread.o lawrite.o  laopen.o laclose.o laio_info.o 
