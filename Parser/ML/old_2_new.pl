#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @a = split(/\t/, $_);

	my $len = length($a[9]);
	my $qual = "";
	
	for (my $i = 0; $i < $len; $i++) {
		$qual .= "I";
	}

	print $a[0], "\t", $a[1], "\t",$a[2], "\t",$a[3], "\t",$a[4], "\t",$a[5], "\t",$a[6], "\t",$a[7], "\t",$a[8], "\t0\t",$a[9], "\t", $qual, "\t", $a[11], "\n";


}

close IN;

exit(0);
