#! /usr/bin/perl
use warnings;
use strict;

my $file = shift;
open FILE, $file or die "cannnot open $file\n";

while(my $line = <FILE>) {
	$line =~ s/^X/23/g;
	$line =~ s/^Y/24/g;
	$line =~ s/^M/25/g;
	print $line;
}

exit(0);
