#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char* argv[]) {
  int i, me, np, nameSize;
  char myProcName[MPI_MAX_PROCESSOR_NAME];
  MPI_Init( &argc , &argv );
  MPI_Comm_rank( MPI_COMM_WORLD , &me );
  MPI_Comm_size( MPI_COMM_WORLD , &np );
  MPI_Get_processor_name( myProcName , &nameSize );
  for( i=0 ; i<np ; ++i ) {
    if( i == me) {
      printf(" rank = %d , processor = %s\n" , me , myProcName);
      fflush(stdout);
    }
    MPI_Barrier( MPI_COMM_WORLD );
  }
  MPI_Finalize();
  exit(0);
}
