#! /usr/bin/perl

use strict;

# The transformation of Col-0 alignments was wrong:
# Just copy the Strain Col-0 alignment to the Reference alignment if different

my $usage = "$0 map.list\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $eco = $a[0];
	if ($eco eq "Col-0") {
		if ($a[4] eq $a[5]) {
			print $line;
		}
		else {
			print $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3], "\t", $a[5], "\t", $a[5];
			for (my $i = 6; $i < @a+0; $i++) {
				print "\t", $a[$i];
			}
			print "\n";
		}
	}
	else {
		print $line;
	}
}

