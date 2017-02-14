#!/usr/bin/perl
use strict;
use Cwd;
use Time::HiRes qw(usleep);   # HACK else > 10 simultaneous ssh's causes some to fail

my($max, %NODES, %PIDS);

unless($max = $ARGV[0])
{
    die "Usage: $0 MAX_SIMULTANEOUS_APPS < APP_LIST\n";
}
if(-t STDIN)    # input not redirected
{
    die "Usage: $0 MAX_SIMULTANEOUS_APPS < APP_LIST\n";
}
print STDERR "Max simultaneous apps: $max\n";

die "Unable to open file PBS_NODEFILE ($ENV{'PBS_NODEFILE'})\n" unless(open(FD, $ENV{'PBS_NODEFILE'}));
while(my $node = <FD>)
{
    chop($node);
    $NODES{$node} = 0;
}
close(FD);
print STDERR "Number of nodes: " . scalar(keys(%NODES)) . "\n";

my $pwd = getcwd;

for(my $i = 0; $i < $max; $i++)
{
    if(my $prog = <STDIN>)
    {
        chop($prog);
        my $node = (sort {$NODES{$a} <=> $NODES{$b}} keys %NODES)[0];
        die "Unable to fork\n" unless(defined(my $pid = fork));
        if($pid == 0)    # child
        {
            die "Unable to execute $prog\n" unless(exec("/usr/bin/ssh -q -n $node \"cd $pwd; $prog\""));
            exit(0);
        }
        $NODES{$node}++;
        $PIDS{$pid} = $node;
#        print STDERR "$prog is $pid on $node\n";
        usleep(100000);
    }
}

while(my $prog = <STDIN>)
{
    chop($prog);
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

for(my $i = 0; $i < $max; $i++)
{
    wait;
}

print STDERR "DONE\n";
