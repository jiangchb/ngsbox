#!/usr/bin/perl
use strict;
use warnings;

my $file    = shift;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	chomp;
	my @e = split(/\t/, $_);

	if(@e != 16) {
		print STDERR "Not enough columns: $e[0]\t$e[1]\n";
	}

	my $brackets = 0;
	for(my $i = 0; $i < length($e[4]); $i++) {

		my $base = substr($e[4], $i, 1);
		if( $base eq "(" || $base eq "[" ) {
			$brackets++;
		}
	}

	print "$e[0]\t$e[1]\t$e[2]\t$e[3]\t$e[4]\t$e[5]\t$e[6]\t$e[7]\t$brackets\t$e[9]\t$e[10]\t$e[11]\t$e[12]\t$e[13]\t$e[14]\t$e[15]\n";
}

close IN;

exit(0);
