#!/usr/bin/perl

use strict;
use warnings;

use threads;
use threads::shared;
use IPC::Open2;
use Thread::Queue;
use List::Util qw(sum);
use File::Basename;
use File::Copy;
use Time::HiRes qw(time);

my $filename = $ARGV[1];
my $nThreads = $ARGV[0];


my $base_filename = basename($filename);
copy($filename,$base_filename); 
$filename = $base_filename;
    
my $process_q = Thread::Queue -> new(); 

my $searchmc = "./SearchMC.pl";

open my $handle, '<', $filename;
chomp(my @lines = <$handle>);
close $handle;

my @targets;
my @upperbounds :shared;
my $start :shared;
my $end :shared;

for my $line (@lines) {
    if ($line =~ /\s*\(declare-fun\s+influence-target-([0-9]*)\s*/) {
        push @targets, "influence-target-$1";
        my $inputfile = "slice-influence-target-$1-$filename";
        $process_q -> enqueue ($inputfile);
    }
}
$process_q -> end();

for my $target_var (@targets) {
    system("./smtslice.pl $target_var $filename");
}
$start = time();
my @jobs = initThreads();

for (my $i = 0; $i < $nThreads; $i++) {
    $jobs[$i] = threads->new(\&runSearchMC);
}

for (my $i = 0; $i < $nThreads; $i++) {
    $jobs[$i]->join();
}
$end = time();
my $upper_bound = sum(@upperbounds);
printf ("Total Sound Upper Bound: %.4f\n", $upper_bound);
printf ("Total Running Time: %.4f s\n", $end - $start);

#print "\n";
#print scalar @upperbounds;
#print "\n";

sub initThreads {
    my @initThreads;
    for(my $i = 1; $i <= $nThreads; $i++) {
        push (@initThreads, $i);
    }
    return @initThreads;
}

sub runSearchMC {
    while (my $filename = $process_q ->dequeue()) {
        my @info = split /-/, $filename;
	printf ("Running command: $searchmc -cl=0.9 -thres=2 -verbose=0 -input_type=smt -solver=cryptominisat -term_cond=1 -output_name=influence-target-$info[3] $filename | grep -v Result\n");
        my $cmd_pid = open2(*OUT, *IN, "$searchmc -cl=0.9 -thres=2 -verbose=0 -input_type=smt -solver=cryptominisat -term_cond=1 -output_name=influence-target-$info[3] $filename | grep -v Result");
        my $line = <OUT>;
        my @result = split ' ', $line;
        close IN;
        close OUT;
        waitpid($cmd_pid, 0);    

        push @upperbounds, $result[2];
	printf ("$filename: Upper Bound is %.4f\n", $result[2]);
        unlink $filename;
        #print scalar @upperbounds;
    }
}
