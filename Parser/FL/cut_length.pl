#!/usr/bin/perl

use strict;
use warnings;

my $usage = "$0 length file\n";
my $length = shift or die $usage;
my $file   = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	
	if( length($elements[1]) >= $length ) {
		my $sequence = substr($elements[1], 0, $length);
		my $prb      = substr($elements[3], 0, $length);
		my $qCal     = substr($elements[4], 0, $length);
		my $chas     = substr($elements[5], 0, $length);
		
		print "$elements[0]\t$sequence\t$elements[2]\t$prb\t$qCal\t$chas\n";
	}
}

close IN;

exit(0);
