#! /usr/bin/perl

use strict;
my $usage = "$0 fasta\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

my $header = "";
while (my $line = <FILE>) {
	chomp($line);
	if (substr($line, 0, 1) ne ">") {
		for (my $i = 0; $i < length($line); $i++) {
			if (	substr($line, $i, 1) ne "A" and substr($line, $i, 1) ne "C" and substr($line, $i, 1) ne "G" and substr($line, $i, 1) ne "T" and
				substr($line, $i, 1) ne "a" and substr($line, $i, 1) ne "c" and substr($line, $i, 1) ne "g" and substr($line, $i, 1) ne "t" and
				substr($line, $i, 1) ne "N"
			) {
				print "$header\n";
				print substr($line, $i, 1) . "\n";
			}
		}
	}
	else {
		$header = $line;
	}
}


