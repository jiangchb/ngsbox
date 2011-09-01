#!/usr/bin/perl

use strict;
use warnings;

my $usage = "$0 file\n";
my $file   = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @e = split("\t", $_);

	if( length($e[1]) != length($e[3]) ) {
		print "$_\n";
	}
}

close IN;

exit(0);
