#! /usr/bin/perl
use strict;

my $usage = "$0 duplication.txt\n";

open FILE, shift or die $usage;

my $sample = "";
my $chr = -1;
my $start = -1;
my $end = -1;
my $cnv_count = 0;

while (<FILE>) {
	my @a = split " ";
	if ($chr != $a[1] or $a[2] > $end + 5000) {
		if ($chr != -1) {
			print $sample, "\t", $chr, "\t", 	$start, "\t", $end, "\t", $cnv_count, "\n";
		}
		$cnv_count = 0;
		$start = $a[2];
	}

	$sample = $a[0];
	$chr = $a[1];
	$end = $a[3];
	$cnv_count += $a[5];
}

if ($chr != -1) {
	print $sample, "\t", $chr, "\t",  $start, "\t", $end, "\t", $cnv_count, "\n";
}



