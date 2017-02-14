#!/bin/ksh

#cd $PBS_O_WORKDIR
cd $JOBDIR

JID=$1
TBEGIN=`echo "print time();" | perl`
./serial_test.x  inp1   > job-stdout0.o$JID
TEND=`echo "print time();" | perl`
echo "Elapsed serial time for 1 instance: `expr $TEND - $TBEGIN`"

TBEGIN=`echo "print time();" | perl`
(   ./serial_test.x  inp1   > job-stdout1.o$JID    ) &
(   ./serial_test.x  inp2   > job-stdout2.o$JID    ) &
(   ./serial_test.x  inp3   > job-stdout3.o$JID    ) &
(   ./serial_test.x  inp4   > job-stdout4.o$JID    ) &
(   ./serial_test.x  inp5   > job-stdout5.o$JID    ) &
(   ./serial_test.x  inp6   > job-stdout6.o$JID    ) &
(   ./serial_test.x  inp7   > job-stdout7.o$JID    ) &
(   ./serial_test.x  inp8   > job-stdout8.o$JID    ) &

###   to prove they are running and on different cores  (PSR)
sleep 5 
echo "executing on `hostname`" 
echo "w command gives"
w
echo ""
echo "------------------------------------------------------------------------"
( ps -e -o user,pid,ppid,lwp,pcpu,psr,etime,vsz,rss,stat,cmd | grep $USER ) || echo "Unable to run ps command"
echo "------------------------------------------------------------------------"

###   wait for all background jobs to end, to show ending times 
wait  

TEND=`echo "print time();" | perl`
echo "Elapsed serial time for 8 instances: `expr $TEND - $TBEGIN`"

echo "  Serial  ended  `date` "
exit
