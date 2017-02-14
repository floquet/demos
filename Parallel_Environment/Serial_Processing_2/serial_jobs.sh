#!/bin/ksh

cd $JOBDIR

TBEGIN=`echo "print time();" | perl`
./serial_test.x  inp1   > job-stdout0.o$1
TEND=`echo "print time();" | perl`
echo "Elapsed serial time for 1 instance: `expr $TEND - $TBEGIN`"

TBEGIN=`echo "print time();" | perl`
(   ./serial_test.x  inp1   > job-stdout1.o$1    ) &
(   ./serial_test.x  inp2   > job-stdout2.o$1    ) &
(   ./serial_test.x  inp3   > job-stdout3.o$1    ) &
(   ./serial_test.x  inp4   > job-stdout4.o$1    ) &
(   ./serial_test.x  inp5   > job-stdout5.o$1    ) &
(   ./serial_test.x  inp6   > job-stdout6.o$1    ) &
(   ./serial_test.x  inp7   > job-stdout7.o$1    ) &
(   ./serial_test.x  inp8   > job-stdout8.o$1    ) &

###   To prove they are running and doing so on different cores (PSR).
###   This step is NOT NECESSARY in your research.

sleep 2 
hostname
w
ps -o user,pid,ppid,lwp,pcpu,psr,etime,vsz,rss,stat,cmd -e | grep -e PID -e ${USER}

###   Wait for all background jobs to end. NECESSARY
wait  

TEND=`echo "print time();" | perl`
echo "Elapsed serial time for 8 instances: `expr $TEND - $TBEGIN`"

exit

