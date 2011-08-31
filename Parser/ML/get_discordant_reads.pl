#! /usr/bin/perl

use strict;

my $usage = "perl map.list";
my $file = shift or die $usage;

open FILE, $file;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[9] == 4 || $a[9] == 5 || $a[9] == 7 || $a[9] == 8 || $a[9] == 10 || $a[9] == 11 || $a[9] == 13 || $a[9] == 14) {
		print $line;
	}
}


