
#ifdef K_AND_R_C

#define name2(a,b) gEnErIc2(a,b)
#define gEnErIc2(a,b) a/**/b

#define name3(a,b,c) gEnErIc3(a,b,c)
#define gEnErIc3(a,b,c) a/**/b/**/c

#define name4(a,b,c,d) gEnErIc4(a,b,c,d)
#define gEnErIc4(a,b,c,d) a/**/b/**/c/**/d

#else /* ANSI_C */

#define name2(a,b) gEnErIc2(a,b)
#define gEnErIc2(a,b) a ## b

#define name3(a,b,c) gEnErIc3(a,b,c)
#define gEnErIc3(a,b,c) a ## b ## c

#define name4(a,b,c,d) gEnErIc4(a,b,c,d)
#define gEnErIc4(a,b,c,d) a ## b ## c ## d

#endif

#ifdef FTN_UNDERSCORES

#define FTN(fname) name2(fname,_)

#else

#define FTN(fname) fname

#endif

