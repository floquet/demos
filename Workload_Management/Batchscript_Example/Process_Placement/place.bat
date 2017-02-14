#!/bin/bash
#PBS -l select=2:ncpus=36:mpiprocs=36
#PBS -l walltime=00:05:00
#PBS -A project-id
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

#  Default placement of processes has first n/p placed on first node,
#  second n/p placed on second node, and so on:
mpiexec_mpt -np 72 ./show_placement.exe

#  However, by obtaining the node names and writing a new PBS_NODEFILE,
#  it is possible to change the process placement.
#  For example, here, we write a node file for a round-robin type of process
#  placement:  
echo "before reset, PBS_NODEFILE is $PBS_NODEFILE"
hosts="`cat $PBS_NODEFILE | uniq`"
let "num_nodes = `cat $PBS_NODEFILE | uniq | wc -l`"
echo "num_nodes = $num_nodes"
let "num_entries = `cat $PBS_NODEFILE | wc -l`" 
let "ii = 0"
touch hostfile.o$JID
until [ "$ii" == "$num_entries" ] ; do
  echo "$hosts" >> hostfile.o$JID
  let "ii = $ii + $num_nodes"  
done
echo "hostfile.o$JID has `cat hostfile.o$JID | wc -l` lines"
TEMPFILE="$PBS_NODEFILE"
export PBS_NODEFILE="`pwd`/hostfile.o$JID"
echo "after PBS_NODEFILE is $PBS_NODEFILE"

mpiexec_mpt -np 72 ./show_placement.exe

#  In order to use Intel MPI's mpirun, swap modules and then use the standard
#  mpirun command:
#  . $MODULESHOME/init/bash
#  module swap mpi/sgimpt mpi/intelmpi
#  module list
#
#  mpirun -np 72 ./show_placement.exe_intel
#
#  The Hydra option -rr under Intel IMPI does not appear to work:
#  mpirun -rr -np 72 ./show_placement.exe
#  nor does -perhost <n>, etc etc

#  Export old value of PBS_NODEFILE in case any post processing system software 
#  depends upon it...
PBS_NODEFILE="$TEMPFILE"

exit
