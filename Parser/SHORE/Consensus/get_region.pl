#! /usr/bin/perl
use strict;

my $usage= "\n$0 consensus_summary.txt chr begin end\n\n" ;
my $file  = shift or die $usage;
my $chr   = shift or die $usage;
my $begin = shift or die $usage;
my $end   = shift or die $usage;

open FILE, $file or die $usage;
my $flag = 0;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[0] == $chr and $a[1] >= $begin and $a[1] <= $end) {
		print $line;
		$flag = 1;
	}
	else {
		if ($flag == 1) {
			exit(1);
		}
	}
}

