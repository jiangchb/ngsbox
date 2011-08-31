#!/usr/bin/perl
use strict;
use warnings;

my $file = shift;
open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;
	my @elem = split(/\t/, $_);

	$elem[10] = "";
	$elem[11] = $elem[9];
	$elem[9]  = 0;


	for (my $i = 0; $i < length($elem[11]); $i++) {
		$elem[10] .= chr(int(64 + 10 * log(1 + 10**((ord(substr($elem[11], $i, 1)) - 33)/10.0)) / log(10)+0.499 ));
		$elem[12] .= "Z";
	}

	print "$elem[0]\t$elem[1]\t$elem[2]\t$elem[3]\t$elem[4]\t$elem[5]\t$elem[6]\t$elem[7]\t$elem[8]\t$elem[9]\t$elem[10]\t$elem[11]\t$elem[12]\n";
}

close IN;

exit(0);
