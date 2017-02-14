#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

int main(int argc, char * argv[]) {
  FILE * pipe;
  char hostname[256];
  pipe = popen("hostname", "r");
  bzero(hostname,256);
  fread(hostname,255,1,pipe);
  pclose(pipe);
  fprintf(stdout,"I am program %s with argument %s on node %s\n",
          argv[0], argv[1], hostname);
  return 0;
}

