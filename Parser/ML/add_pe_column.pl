#!/usr/bin/perl
use strict;
use warnings;

my $file    = shift;
open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	if(@elements == 12) {
		print "$elements[0]\t$elements[1]\t$elements[2]\t$elements[3]\t$elements[4]\t$elements[5]\t$elements[6]\t$elements[7]\t$elements[8]\t0\t$elements[9]\t$elements[10]\t$elements[11]\n";
	}
	else {
		print "$_\n";
	}
}

close IN;

exit(0);
