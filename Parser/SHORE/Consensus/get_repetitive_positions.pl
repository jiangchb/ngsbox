#! /usr/bin/perl
use strict;

my $usage= "\n$0 consensus\n\n" ;
my $file = shift or die $usage;

open FILE, $file or die $usage;

while (<FILE>) {
	chomp;
	my @a = split("\t", $_);

	if( $a[10] > 1 ) {
		print $a[1] . "\t" . $a[10] . "\n";
	}
}
