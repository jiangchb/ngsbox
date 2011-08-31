#! /usr/bin/perl
use strict;

my $usage = "$0 file chr,pos-column chr start end\n";
my $file = shift or die $usage;
my $chrpos = shift or die $usage;
my $chr = shift or die $usage;
my $start = shift or die $usage;
my $end = shift or die $usage;

my ($chrcolumn, $poscolumn) = split ",", $chrpos;
$chrcolumn--;
$poscolumn--;

open FILE, $file or die "Cannot open file\n";

my $flag = 0;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[$chrcolumn]==$chr and $a[$poscolumn]>=$start and $a[$poscolumn]<=$end) {
		$flag = 1;
		print $line;
	}
	else {
		if ($flag == 1) {
			exit(0);
		}
	}
}


