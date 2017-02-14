
#ifndef _MEASURE_MEMORY_HH
#define _MEASURE_MEMORY_HH

/** @file
    @brief A Documented file.
    
    Details.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>

#if UNDERSCORECOUNT == 1
#  define MEASURE_MEMORY measure_memory_
#  define PID_VM pid_vm_
#  define PID_RSM pid_rsm_
#elif UNDERSCORECOUNT == 2
#  define MEASURE_MEMORY measure_memory__
#  define PID_VM pid_vm__
#  define PID_RSM pid_rsm__
#else
#  define MEASURE_MEMORY measure_memory
#  define PID_VM pid_vm
#  define PID_RSM pid_rsm
#endif

#if defined (__cplusplus) || defined (c_plusplus)
extern "C" {
#endif

int MEASURE_MEMORY(const char * key, int64_t *amount);
int PID_VM(const pid_t * pid, int64_t * amount);
int PID_RSM(const pid_t * pid, int64_t * amount);

#if defined (__cplusplus) || defined (c_plusplus)
}
#endif

#endif // _MEASURE_MEMORY_HH

