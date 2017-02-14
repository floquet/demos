#ifndef BITVEC_H
#define BITVEC_H 1

typedef unsigned int bitvector;

extern bitvector *newbitvec( int nbits );

extern void freebitvec( bitvector *bitvec );

extern void setbit( int ithbit, bitvector bitvec[] );


extern void clrbit( int ithbit, bitvector bitvec[] );

extern int testbit( int ithbit, bitvector bitvec[] );


#endif
