#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);
use Data::Dumper qw(Dumper);

my $filename = $ARGV[1];
my $target = $ARGV[0];
my @vars;

open my $handle, '<', $filename;
chomp(my @lines = <$handle>);
close $handle;

#init
my $index=0;
my $line;
my @taint_lines;

for $line (@lines) {
    $taint_lines[$index] = 0;
    if ($line =~ /\s*\(declare-fun\s+(\S+)\s*/) {
        my $var = $1;
        push @vars, $var;
    } elsif ($line =~ /\s*\(assert\s+/){
        if ($line =~ /\s*\(assert\s+true/){
            $taint_lines[$index] = 1;
        } elsif ($line =~ /\s*\(assert\s+true/){
            $taint_lines[$index] = 1;
        }
    } else {
        $taint_lines[$index] = 1;
    }
    $index++;
}

my @taint_vars;
push @taint_vars, $target;
my @new_taint_vars;
my @cond_taint_vars;

#loop for finding tainted lines
FOO1: {
    while () {
        my $a = sum(@taint_lines);
        $index=0;
        for $line (@lines) {
            if ($taint_lines[$index] == 0) {
                for my $taint_var (@taint_vars) {
                    if ($line =~ /\s$taint_var\s/) {
                        $taint_lines[$index] = 1;
                        push @new_taint_vars, parse_vars($line);
                        push @cond_taint_vars, parse_cond_vars($line);
                        @new_taint_vars = uniq(@new_taint_vars);
                        @cond_taint_vars = uniq(@cond_taint_vars);
                    }
                }
            }
            $index++;
        }
        @taint_vars = @new_taint_vars;
        #print join(", ", @new_taint_vars);
        #print "\n";
        my $b = sum(@taint_lines);
        if ($a == $b) {
            last FOO1;
        }
    }
}

#print join(", ", @taint_vars);
#print "\n";
#print join(", ", @cond_taint_vars);
#print "\n";

#condition variables
FOO2: {
    while () {
        my $a = sum(@taint_lines);
        $index=0;
        for $line (@lines) {
            if ($taint_lines[$index] == 0) {
                for my $cond_taint_var (@cond_taint_vars) {
                    if ($line =~ /ite\s$cond_taint_var\s/) {
                        
                    } elsif ($line =~ /\s$cond_taint_var\s/) {
                        $taint_lines[$index] = 1;
                        push @taint_vars, parse_vars($line);
                        push @cond_taint_vars, parse_cond_vars($line);
                        @taint_vars = uniq(@taint_vars);
                        @cond_taint_vars = uniq(@cond_taint_vars);
                    }
                }
            }
            $index++;
        }
        my $b = sum(@taint_lines);
        if ($a == $b) {
            last FOO2;
        }
    }
}

#print join(", ", @taint_vars);
#print "\n";
#print join(", ", @cond_taint_vars);
#print "\n";

FOO3: {
    while () {
        my $a = sum(@taint_lines);
        $index=0;
        for $line (@lines) {
            if ($taint_lines[$index] == 0) {
                for my $taint_var (@taint_vars) {
                    if ($line =~ /\s$taint_var\s/) {
                        $taint_lines[$index] = 1;
                        push @new_taint_vars, parse_vars($line);
                        push @cond_taint_vars, parse_cond_vars($line);
                        @new_taint_vars = uniq(@new_taint_vars);
                        @cond_taint_vars = uniq(@cond_taint_vars);
                    }
                }
            }
            $index++;
        }
        @taint_vars = @new_taint_vars;
        #        print join(", ", @new_taint_vars);
        #        print "\n";
        my $b = sum(@taint_lines);
        if ($a == $b) {
            last FOO3;
        }
    }
}

#print join(", ", @taint_vars);
#print "\n";
#print join(", ", @cond_taint_vars);
#print "\n";

sub parse_vars {
    my ($input_line) = @_;
    my @parsed_vars;
    my @words = split / |\(|\)/, $input_line;
    for my $word (@words) {
        if ($word =~ /^cond/) {
        } elsif ($word =~ /^[-A-Za-z]/) {
            if ( grep( /^$word$/, @vars ) ) {
                push @parsed_vars, $word;
            }
        }
    }
    return @parsed_vars;
}

sub parse_cond_vars {
    my ($input_line) = @_;
    my @parsed_vars;
    my @words = split / |\(|\)/, $input_line;
    for my $word (@words) {
        if ($word =~ /^cond/) {
            if ( grep( /^$word$/, @vars ) ) {
                push @parsed_vars, $word;
            }
        }
    }
    return @parsed_vars;
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

$index=0;
open(my $fh, '>', "slice-$target-$filename");
for $line (@lines) {
    if ($taint_lines[$index] == 1) {
        print $fh "$line\n";
        #print "$line\n";
    }
    $index++;
}
#    for $taint_var (@new_taint_vars) {
#        if ($line =~ /$taint_var/) {
#            print $fh "$line\n";
#            print "$line\n";
#            last;
#        }
#    }
close $fh;

