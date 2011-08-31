#! /usr/bin/perl
use strict;
my $usage = "$0 inverionsprediction\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;
my $c = 0;
my $e = "";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] ne $e) {
		$c++;
		$e = $a[0];
	}
	print $c, "\t", $e, "\t", $a[2], "\t", $a[4], "\t", $a[5], "\t", $a[5]-$a[4]+1, "\n"; 
}
