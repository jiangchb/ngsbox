#!/usr/bin/perl

use strict;
use warnings;

my $file   = shift or die "$0 file";
my $current_lane = -1;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);

	if($current_lane != $e[2]) {
		if($current_lane != -1) { close OUT; }
		$current_lane = $e[2];
		open OUT, ">$file.$current_lane" or die "Cannot open output file";
	}
	print OUT "$_\n";
}

close IN;
close OUT;

exit(0);
