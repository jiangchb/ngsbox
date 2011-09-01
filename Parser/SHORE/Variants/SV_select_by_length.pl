#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 SVdeletion minlength\n\n";
my $file = shift or die $usage;
my $length = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[6] >= $length) {
		print $line;
	}
}



