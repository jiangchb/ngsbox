#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;

my $qual_offset = 0;
my $qual2sanger = 1;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my ($id, $sequence, $pe, $qual1, $qual2, $qual3) = split(/\t/, $_);

	if ($qual_offset == 1) {
		my $new = "";
		for (my $i = 0; $i<length($qual1); $i++) {
			$new .= chr((ord(substr($qual1, $i, 1))+14));
		}
		$qual1 = $new;
	}

	if ($qual2sanger == 1) {
		$qual2 = "";
		for (my $i = 0; $i<length($qual1); $i++) {
			$qual2 .= chr(int(33 + 10 * log(1 + 10**((ord(substr($qual1, $i, 1)) - 64)/10.0)) / log(10)+0.499 ));
		}
	}

	print "$id\t$sequence\t$pe\t$qual1\t$qual2\t$qual3\n";
}

close IN;

exit(0);
