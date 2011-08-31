#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\n$0 maxhit file\n\n";
my $max_hit = shift or die $usage;
my $file    = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	my @elements = split(/\t/, $_);

	if(@elements < 11) {
		print STDERR "Not enough columns: $elements[0]\t$elements[1]\n";
	}
	elsif( $elements[6] !~ /\d+/ ) {
		print STDERR "Hits not int: $elements[0]\t$elements[1]\t$elements[6]\n";
	}
	elsif($elements[6] <= $max_hit) {
		print $_;
	}
}

close IN;

exit(0);
