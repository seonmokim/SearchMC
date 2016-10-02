#!/usr/bin/perl

use strict;
use warnings;

# sub fact {
#     my($n) = @_;
#     my $r = 1;
#     for my $x (2 .. $n) {
# 	$r *= $x;
#     }
#     return $r;
# }

# sub choose {
#     my($n, $k) = @_;
#     return fact($n) / ( fact($n - $k) * fact($k) );
# }

my @pascal;
for my $n (0 .. 1100) {
    $pascal[$n][0] = $pascal[$n][$n] = 1;
    for my $k (1 .. $n - 1) {
	$pascal[$n][$k] = $pascal[$n-1][$k-1] + $pascal[$n-1][$k];
    }
}

sub choose {
    my($n, $k) = @_;
    die "Out of range in choose($n, $k)" unless defined $pascal[$n][$k];
    return $pascal[$n][$k];
}

# for my $n (0 .. 1000) {
#     for my $k (0 .. $n) {
# 	print choose($n, $k), " ";
#     }
#     print "\n";
# }

# exit;

#for my $s (1 .. 16) {
#    my $b = log($s)/log(2);
for my $bt (0, 10, 16, 20, 23, 26, 28, 30, 32, 33, 34, 36 .. 100) {
    my $b = $bt / 10;
    my $s = int(2**$b + 0.5);
    printf "%g bits, $s solutions\n", $b;
    for my $k (0 .. 20) {
	printf "%2d", $k;
	my $rest = 1;
	for my $n (0 .. 15) {
	    last if $n > $s;
	    my $p1 = 2**-$k;
	    my $p2 = 1-$p1;
	    my $p = choose($s, $n) * ($p1 ** $n) * ($p2 ** ($s - $n));
	    $rest -= $p;
	    printf " %.4f", $p;
	}
	printf " %.4f", abs($rest);
	print "\n";
	
    }
    print "\n";
}

