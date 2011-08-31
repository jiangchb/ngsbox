#!/usr/bin/perl
use strict;
use warnings;

my $max_mm = shift;
my $file    = shift;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	my @elements = split(/\t/, $_);

	if(@elements != 13) {
		print STDERR "Not enough columns: $elements[0]\t$elements[1]\n";
	}
	elsif( $elements[5] !~ /\d+/ ) {
		print STDERR "MM not int: $elements[0]\t$elements[1]\t$elements[6]\n";
	}
	elsif($elements[5] <= $max_mm) {
		print $_;
	}
}

close IN;

exit(0);
