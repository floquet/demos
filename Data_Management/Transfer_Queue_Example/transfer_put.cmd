#!/bin/bash
##Here is a sample standalone command file that you may use as a template:
#PBS -l walltime=00:##:00
#PBS -l select=1:ncpus=1
#PBS -A ABCDE00000000
#PBS -q transfer
#PBS -N transfer_put
#PBS -j oe

## See ERDC DSRC Archive User Guide at
## http://centers.hpc.mil/users/documentation.html 
## regarding ERDC DSRC archive policies.

## Pathname to directory to hold the data.
## Edit to your needs
TRANSFER_DIR=$ARCHIVE_HOME/transfer/$PBS_JOBID

## $TRANSFER_DIR not assumed to exist prior to this job, so it must be created
## on the archive system.

archive mkdir -p $TRANSFER_DIR

echo "starting transfer"

cd $WORKDIR

## Please adjust your tar command accordingly, you should use a multipart
## tar file if it exceeds 500GB

tar -pzcvf file_name.tgz <dir_name/file_name>

archive put -C $TRANSFER_DIR file_name.tgz

exit
