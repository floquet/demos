#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
int main(int argc, char* argv[]) {
  int i, me, you_plus, you_minus, np, nameSize;
  char myProcName[MPI_MAX_PROCESSOR_NAME];
  int sendbuf[3], recvbuf[3];
  int tag = 2;
  MPI_Status mpistat;
  MPI_Init( &argc , &argv );
  MPI_Comm_rank( MPI_COMM_WORLD , &me );
  MPI_Comm_size( MPI_COMM_WORLD , &np );
  MPI_Get_processor_name( myProcName , &nameSize );
  for( i=0 ; i<np ; ++i ) {
    if( i == me) {
      printf("part_three.exe: rank = %d , processor = %s\n" , me , myProcName);
      fflush(stdout);
    }
    MPI_Barrier( MPI_COMM_WORLD );
  }
  sendbuf[0] = me;
  recvbuf[0] = -1;
  you_plus = me + 1;
  if(you_plus >= np) { you_plus -= np; }
  you_minus = me - 1;
  if(you_minus < 0) { you_minus += np; }
  if(me != 0) {
    MPI_Send(sendbuf, 1, MPI_INT, you_plus, tag, MPI_COMM_WORLD);
    MPI_Recv(recvbuf, 1, MPI_INT, you_minus, tag, MPI_COMM_WORLD, &mpistat);
  }
  else {
    MPI_Recv(recvbuf, 1, MPI_INT, you_minus, tag, MPI_COMM_WORLD, &mpistat);
    MPI_Send(sendbuf, 1, MPI_INT, you_plus, tag, MPI_COMM_WORLD);
  }
  if(recvbuf[0] != you_minus) {
    fprintf(stderr,"Error: recvbuf[0] = %d , you_minus = %d\n", recvbuf[0], you_minus);
  }
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
  return 0;
}
