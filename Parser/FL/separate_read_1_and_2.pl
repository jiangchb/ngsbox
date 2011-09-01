#!/usr/bin/perl

use strict;
use warnings;

my $file   = shift;

open IN, $file or die "Cannot open input file\n";
open OUT1, ">$file.r1" or die "Cannot open output file\n";
open OUT2, ">$file.r2" or die "Cannot open output file\n";

while( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);
	
	if($e[2] == 1) {
		print OUT1 "$_\n";
	}
	elsif($e[2] == 2) {
		print OUT2 "$_\n";
	}
}

close IN;

exit(0);
