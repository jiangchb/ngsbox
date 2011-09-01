#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\n$0 min_size segment_file\n\n";
my $min_size = shift or die $usage;
my $segment = shift or die $usage;

open SEG, $segment or die "Cannot open input file\n";

while( my $line = <SEG> ) {
	chomp($line);

	if($line =~ /^#/) {
		#print "$line\n";
	}
	else {
		my @a = split("\t", $line);
		my $beg = $a[1] - 1;
		my $end = $beg + $a[2];
		if($a[2] >= $min_size) {
			print "$a[0]\t$beg\t$end\n";
		}
	}
}

close SEG;

exit(0);
