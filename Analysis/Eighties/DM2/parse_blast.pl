#! /usr/bin/perl

use strict;
use warnings;

my $usage = "$0 blnfile\n";

my $file = shift or die $usage;
open FILE, $file or die $usage;

my $query = "";
my $max = -1;
my $min = -1;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (substr($line, 0, 6) eq "Query=") {
		$query = $a[1];
	}
	elsif (substr($line, 0, 8) eq " Score =" || substr($line, 0, 14) eq "Gap Penalties:") {
		if ($min != -1) {
			if ($max - $min > 300) {
				print $query, "\t", $min, "\t", $max, "\n";
			}
			$min = -1;
			$max = -1;
		}
	}
	elsif (substr($line, 0, 6) eq "Query:") {
		if ($min == -1 || $a[1] < $min) {
			$min = $a[1];
		}
		if ($max == -1 || $a[1] > $max) {
			$max = $a[1];
		}
		if ($min == -1 || $a[3] < $min) {
                        $min = $a[3];
                }
                if ($max == -1 || $a[3] > $max) {
                        $max = $a[3];
                }
	}
}


