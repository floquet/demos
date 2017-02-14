#!/bin/bash
#PBS -l select=2:ncpus=36:mpiprocs=36
#PBS -l place=scatter:excl:group=switch
#PBS -l walltime=00:05:00
#PBS -A <your-project-id>
#PBS -q debug
#PBS -N place
#PBS -j oe
#PBS -l application=other

cd $PBS_O_WORKDIR

# Set up execution directory in $WORKDIR
JID=`echo $PBS_JOBID | cut -d. -f1`
export JOBDIR=$WORKDIR/placement_o$JID
mkdir $JOBDIR

# Copy binary to execution directory and move into it
cd $PBS_O_WORKDIR
cp show_placement.exe $JOBDIR
cd $JOBDIR

mpiexec_mpt -np 72 ./show_placement.exe

#  In order to use mpirun under SGI MPT, set up the hosts as a comma-separated
#  list, using $PBS_NODEFILE contents, and then use a modified mpirun command.
#  Otherwise, mpirun tries to place all the MPI processes on the first node
#  only.
#  hosts="`cat $PBS_NODEFILE | uniq`"
#  par_host=""
#  for node in $hosts ; do
#    par_host="${par_host}, ${node}"
#  done
#  par_host=`echo $par_host | sed "s/^, //"`
#  
#  The following command will place 36 MPI processes on each host in $par_host:
#  mpirun $par_host -np 36 ./show_placement.exe


#  In order to use Intel MPI's mpirun, swap modules and then use the standard
#  mpirun command.  This requires the binary to be linked against Intel's IMPI.
#  . $MODULESHOME/init/bash
#  module swap mpi/sgimpt mpi/intelmpi
#  module list
#
#  mpirun -np 72 ./show_placement.exe

exit
