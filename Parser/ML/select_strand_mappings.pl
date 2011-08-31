#! /usr/bin/perl


use strict;
my $usage = "perl $0 maplist D|P";

my $file = shift or die $usage;
my $dir = shift or die $usage;
die $usage if ($dir ne "P" and $dir ne "D");

open FILE, $file or die $usage;

while (<FILE>) {
	my @a = split " ";
	if ($a[4] eq $dir) {
		print $_;
	}
}


