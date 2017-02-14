#include <stdio.h>
#include <stdlib.h>
#include <omp.h>           // needed for OpenMP
#include <unistd.h>        // only needed for definition of gethostname
#include <sys/param.h>     // only needed for definition of MAXHOSTNAMELEN

int main (int argc, char *argv[]) {
  int th_id, nthreads;
  char foo[] = "Hello";
  char bar[] = "World";
  char hostname[MAXHOSTNAMELEN];
  gethostname(hostname, MAXHOSTNAMELEN);
  #pragma omp parallel private(th_id)
  {
    th_id = omp_get_thread_num();
    printf("%s %s from thread %d on %s!\n", foo, bar, th_id, hostname);
    #pragma omp barrier
    if ( th_id == 0 ) {
      nthreads = omp_get_num_threads();
      printf("There were %d threads on %s!\n", nthreads, hostname);
    }
  }
  return EXIT_SUCCESS;
}
