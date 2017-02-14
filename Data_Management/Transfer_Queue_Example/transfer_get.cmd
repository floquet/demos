#!/bin/bash
##Here is a sample standalone command file that you may use as a template:
#PBS -l walltime=00:##:00
#PBS -l select=1:ncpus=1
#PBS -A ABCDE00000000
#PBS -N transfer_get
#PBS -q transfer
#PBS -j oe

## See ERDC DSRC Archive User Guide at
## http://centers.hpc.mil/users/documentation.html 
## regarding ERDC DSRC archive policies.

## TRANSFER_DIR holds the data, so exists already!
## Edit to your needs
TRANSFER_DIR=$ARCHIVE_HOME/transfer/$PBS_JOBID

echo "starting transfer"

## Be certain to transfer to temporary workspace to avoid violating quotas
## on $HOME

cd $WORKDIR

## You may need to adjust your archive and tar commands accordingly if
## the archive tar files exceed 500GB

archive get -C $TRANSFER_DIR file_name.tgz

tar -xzvf file_name.tgz 

exit
