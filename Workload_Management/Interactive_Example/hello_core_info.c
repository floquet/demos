#define _GNU_SOURCE
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#include <mpi.h>
#include <omp.h>

/* On SGI systems: */
/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c and Cray's xthi.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return(str);
}

/* On Cray systems: */
/* int get_core_num()
{
    int core;
    cpu_set_t coremask;
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);

    for (core = 0; core < CPU_SETSIZE; core++)
    {
        if (CPU_ISSET(core, &coremask))
        {
          return(core);
        }
    }
    return(-1);
}

int get_core_num_() { return(get_core_num()); } */

void print_hello_info(MPI_Comm comm)
{
    int rank, thread;
    cpu_set_t coremask;
    char core[7 * CPU_SETSIZE];
    char hostname[HOST_NAME_MAX];

    memset(core, 0, sizeof(core));
    (void)MPI_Comm_rank(comm, &rank);
    (void)gethostname(hostname, HOST_NAME_MAX);
    {
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, core);
        /* core = get_core_num(); */
        printf("Hello from rank %2d, host %s, core %s\n", rank, hostname, core);
    }
}

void print_hello_info_(MPI_Comm *comm) { print_hello_info(*comm); return; }

/* End of core_info */
