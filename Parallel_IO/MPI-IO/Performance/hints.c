
/*
 * Get hints from an open file.
 */
/*
 * Code taken from IOR which is copyrighted by
 * The Regents of the University of California
 * Their license restrictions are the following.
 * You can redistribute it and/or modify it under
 * the terms of the GNU General Public License
 * (as published by the Free Software
 * Foundation) version 2, dated June 1991.
 */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#ifndef MPI_FILE_NULL
#   include <mpio.h>
#endif /* not MPI_FILE_NULL */

void ShowHints (MPI_Info *);

/*
 * MPI_CHECK will display a custom error message as well as an error string
 * from the MPI_STATUS and then exit the program
 */

#define MPI_CHECK(MPI_STATUS, MSG) do {                                  \
    char resultString[MPI_MAX_ERROR_STRING];                             \
    int resultLength;                                                    \
                                                                         \
    if (MPI_STATUS != MPI_SUCCESS) {                                     \
        fprintf(stdout, "** error **\n");                                \
        fprintf(stdout, "ERROR in %s (line %d): %s.\n",                  \
                __FILE__, __LINE__, MSG);                                \
        MPI_Error_string(MPI_STATUS, resultString, &resultLength);       \
        fprintf(stdout, "MPI %s\n", resultString);                       \
        fprintf(stdout, "** exiting **\n");                              \
        fflush(stdout);                                                  \
        MPI_Abort(MPI_COMM_WORLD, -1);                                   \
    }                                                                    \
} while(0)


int main(int argc, char * argv[]) {

  int num_tasks, rank;
  MPI_Info mpiHints = MPI_INFO_NULL;
  MPI_File   * fd;

  char testFileName[] = "test_file.tmp";

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  fd = (MPI_File *) malloc(sizeof(MPI_File));

  *fd = 0;
  MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, testFileName,
			  MPI_MODE_RDWR, mpiHints, fd),
	    "cannot open file");

  if (rank == 0 ) {
    MPI_CHECK(MPI_File_get_info(*fd, &mpiHints),
	      "cannot get file info");
    fprintf(stdout, "\nhints passed to MPI_File_open() {\n");
    ShowHints(&mpiHints);
    fprintf(stdout, "}\n");
  }

  MPI_File_close((MPI_File *)fd);
  MPI_Finalize();

  return 0;

} // main


/*
 * Show all hints (key/value pairs) in an MPI_Info object.
 */

void ShowHints(MPI_Info * mpiHints)
{
    char key[MPI_MAX_INFO_VAL],
         value[MPI_MAX_INFO_VAL];
    int  flag,
         i,
         nkeys;

    MPI_Info_get_nkeys(*mpiHints, &nkeys);
    /*
    MPI_CHECK(MPI_Info_get_nkeys(*mpiHints, &nkeys),
              "cannot get info object keys");
    */

    for (i = 0; i < nkeys; i++) {
        MPI_CHECK(MPI_Info_get_nthkey(*mpiHints, i, key),
                  "cannot get info object key");
        MPI_CHECK(MPI_Info_get(*mpiHints, key, MPI_MAX_INFO_VAL-1,
                               value, &flag),
                  "cannot get info object value");
        fprintf(stdout,"\t%s = %s\n", key, value);
    }
} /* ShowHints() */

