#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;
open FILE, $file or die "Cannot open infile\n";

while(<FILE>) {
	my $line = $_;
	my @elem = split("\t", $line);

	if( 	($elem[6] < 0.0000000001) && 
		($elem[5] >= 5) &&
		($elem[4] < 75) &&
		($elem[3] < 0)
	) {
		print "$line";
	}
}

exit(0);
