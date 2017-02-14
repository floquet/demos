#!/usr/bin/perl -Wall
use strict;
use Cwd;
use Time::HiRes qw(usleep);   # HACK else > 10 simultaneous ssh's causes some to fail

my($max, %NODES, %PIDS);

$\ = ""; # No extra carriage return.

# Read run parameters

unless($max = $ARGV[0])
{
    print STDERR "Error, maximum simultaneous applications not specified.\n";
    die "Usage: $0 MAX_SIMULTANEOUS_APPS < APP_LIST\n";
}
if(-t STDIN)    # input not redirected
{
    print STDERR "Error, no stdin file specified.\n";
    die "Usage: $0 MAX_SIMULTANEOUS_APPS < APP_LIST\n";
}
print STDERR "Max simultaneous apps: $max\n";

# Open file specified by hostlist.
# Sometimes, after maintenance, $PBS_NODEFILE is not available
# on the compute nodes.
# On the batch node, which interprets the pack.pbs script,
# copy $PBS_NODEFILE to file hostlist.

# Read list of compute nodes.
# Hash %NODES has key "$node" only once, though it would appear
# more than once in hostlist file.

local *FD;
open(FD,'<','hostlist') || die 'Unable to open hostlist';
my $compute_node;
while(<FD>) {
    $compute_node = $_;
    chomp($compute_node); # Get rid of final carraige return.
    $NODES{$compute_node} = 0;
}
close(FD);
print STDERR "Number of nodes: " . scalar(keys(%NODES)) . "\n";

my $pwd = getcwd;

# Start first $max processes.
for(my $i = 0; $i < $max; $i++)
{
    # Read next item in stdin, if there is one.
    # <STDIN> reads one line.
    # Reading beyond last item returns null.
    my $prog = <STDIN>;
    if(defined($prog) && ($prog ne ""))
    {
        chomp($prog); # Get rid of final carraige return.
	# Sort keys of %NODES, which are node names.
        # Use for sorting the values of $NODES{some_node}.
        # Assign to $node the first of the sorted list, which
        # corresponds to the lowest value.
        # The fact that there are 32 cores per node is not a
        # consideration for this script.  The tasks are distributed
        # uniformly up to $max at any one time.
        my $node = (sort {$NODES{$a} <=> $NODES{$b}} keys %NODES)[0];
        # Fork, both parent and child execute this script.
        die "Unable to fork\n" unless(defined(my $pid = fork));
        if($pid == 0)    # child
        {
            # The chosen $node executes the composite command of
            # changing to the current directory and then executing
            # the text read from one line of STDIN.
            die "Unable to execute $prog\n" unless(exec("/usr/bin/ssh -q -n $node \"cd $pwd; $prog\""));
            exit(0);
        }
        # Only the parent arrives here because after the child executes
        # remote task using ssh it encounters exit(0).
        # (A "remote" task could be on node in which parent is
        # running this script.  This script is being run on a
        # compute node, not on the batch node.)
        # The value in hash $NODE associated with the node $node
        # increases to indicate that $node has been assigned a task.
        $NODES{$node}++;
        # The PID of the child process used as a hash key,
        # with the value being $node.
        $PIDS{$pid} = $node;
#        print STDERR "$prog is $pid on $node\n";
        # A tenth of a second delay between calls to ssh.
        usleep(100000);
    }
}

# Read next item in stdin, while there is one.
# Reading beyond last item returns null.
# Algorithm below is same as above except that
# (1) each "while" cycle starts when a forked process
#     finishes resulting in the wait command returning.
# (2) the count of active processes in %NODES is decremented
#     and the hash item $PIDS{$kid} is removed.
while(my $prog = <STDIN>)
{
    chomp($prog); # Get rid of final carraige return.
    # Wait for a process to finish.
    my $kid = wait;
#    print STDERR "$kid DIED " . ($? >> 8) . "\n";
    $NODES{ $PIDS{$kid} }--;
    delete $PIDS{$kid};

    my $node = (sort {$NODES{$a} <=> $NODES{$b}} keys %NODES)[0];
    die "Unable to fork\n" unless(defined(my $pid = fork));
    if($pid == 0)    # child
    {
        die "Unable to execute $prog\n" unless(exec("/usr/bin/ssh -q -n $node \"cd $pwd; $prog\""));
        exit(0);
    }
    $NODES{$node}++;
    $PIDS{$pid} = $node;
#    print STDERR "$prog is $pid on $node\n";
    usleep(100000);
}

# Wait for all forked processes to finish.
# If parent process exits before child processes finish,
# then the child processes will become zombies.

my $kiddie;
do {
  $kiddie = waitpid(-1, 0);
} while $kiddie > 0;

print STDERR "DONE\n";

exit 0;

__END__


