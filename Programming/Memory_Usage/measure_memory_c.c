
#include "measure_memory.h"

#if defined (__cplusplus) || defined (c_plusplus)
extern "C" {
#endif

int MEASURE_MEMORY(const char * key, int64_t *amount) {
  int status = 0;
  int i;
  char * got;
  int scan_count;
  char oneline[512];
  char * fgets_retval;
  int64_t amount_tmp;
  FILE* vmstat_pipe = popen("/usr/bin/vmstat -s","r");
  if(vmstat_pipe == NULL) {
    fprintf(stderr,"Error in executing /usr/bin/vmstat -s\n");
    fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
  }
  else {
    fgets_retval = oneline;
    while(fgets_retval != NULL) {
      for(i = 0; i < 20; i++) {
	fgets_retval = fgets(oneline, 511, vmstat_pipe);
	if(fgets_retval == NULL) { break; }
	got = strstr(oneline, key);
	if(got != NULL) {
	  scan_count = sscanf(oneline, "%ld", &amount_tmp);
	  if(scan_count == 1) {
	    *amount = amount_tmp;
	  }
	  else {
	    fprintf(stderr,"Error reading free memory integer\n");
	    fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
	  }
	}
      }
    }
    pclose(vmstat_pipe);
  }
  return status;
} // MEASURE_MEMORY

int PID_VM(const pid_t * pid, int64_t * amount) {
  int status = 0;
  char command[512];
  bzero(command, 512);
  int64_t amount_tmp;
  int scan_count;
  FILE* vmstat_pipe;
  sprintf(command, "/bin/ps --pid %d --format vsize --no-headers",
	  *pid);
  vmstat_pipe = popen(command,"r");
  if(vmstat_pipe == NULL) {
    fprintf(stderr,"Error in executing ps --format vsize\n");
    fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
  }
  else {
    scan_count = fscanf(vmstat_pipe, "%ld", &amount_tmp);
    if(scan_count == 1) {
      *amount = amount_tmp;
    }
    else {
      fprintf(stderr,"Error reading pid resident memory integer\n");
      fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
    }
    pclose(vmstat_pipe);
  }
  return status;
} // PID_VM

int PID_RSM(const pid_t * pid, int64_t * amount) {
  int status = 0;
  char command[512];
  bzero(command, 512);
  int64_t amount_tmp;
  int scan_count;
  FILE* rsmstat_pipe;
  sprintf(command, "/bin/ps --pid %d --format rssize --no-headers",
	  *pid);
  rsmstat_pipe = popen(command,"r");
  if(rsmstat_pipe == NULL) {
    fprintf(stderr,"Error in executing ps --format rssize\n");
    fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
  }
  else {
    scan_count = fscanf(rsmstat_pipe, "%ld", &amount_tmp);
    if(scan_count == 1) {
      *amount = amount_tmp;
    }
    else {
      fprintf(stderr,"Error reading pid resident memory integer\n");
      fprintf(stderr,"line %d file %s\n", __LINE__, __FILE__);
    }
    pclose(rsmstat_pipe);
  }
  return status;
} // PID_RSM

#if defined (__cplusplus) || defined (c_plusplus)
}
#endif
