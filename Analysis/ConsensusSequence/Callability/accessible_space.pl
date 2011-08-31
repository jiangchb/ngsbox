#!/usr/bin/perl

use strict;
use warnings;

my $file = shift or die;

my $callable = 0;
my $not_callable = 0;

open FILE, $file or die "Cannot open input file.\n";
while( <FILE> ) {
	chomp;
	my @a = split("\t", $_);
	if($a[5] >= 25 && $a[7] >= 0.7) {
		$callable++;
	}
	else {
		$not_callable++;
	}
}

print "$callable\t$not_callable\n";

exit(0);
