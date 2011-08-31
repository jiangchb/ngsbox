#! /usr/bin/perl
use strict;

my $usage = "$0 file(chr pos)\nWill be translated in (chr begin end)\n";
my $file = shift or die $usage;

open FILE, $file or die "Cannot open file\n";

my $bc = -1;
my $bp = -1;
my $lp = -1;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] > $bc or $a[1] > $lp+1) {
		if ($bc != -1) {
			print $bc, "\t", $bp, "\t", $lp, "\n";
		}
		$bc = $a[0];
		$bp = $a[1];
	}
	$lp = $a[1];
}
if ($bc != -1) {
	print $bc, "\t", $bp, "\t", $lp, "\n";
}

