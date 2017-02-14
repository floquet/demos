
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILE_DESC 16384

int main(int argc, char * argv[]) {
   FILE * fd[MAX_FILE_DESC];
   char filename[MAX_FILE_DESC][512];
   int i;
   for(i = 0; i < MAX_FILE_DESC; i++) {
      fd[i] = NULL;
      sprintf(filename[i],"aaastub%d",i);
   }
   for(i = 0; i < MAX_FILE_DESC; i++) {
      fd[i] = fopen(filename[i],"w");
      if(fd[i] == NULL) {
         fprintf(stdout,"Error, unable to open %s\n",
                 filename[i]);
         return 1;
      }
   }
   for(i = 0; i < MAX_FILE_DESC; i++) {
      if(fd[i] != NULL) { fclose(fd[i]); }
   }


   return 0;
}
