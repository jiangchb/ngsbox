#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\n$0 file\n\n";
my $file  = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	my @e = split(/\t/, $_);

	if( $e[2] =~ /\[L/ ) {
		$e[1]++;
		print $e[0]."\t".$e[1]."\t".$e[2]."\t".$e[3]."\t".$e[4]."\t".$e[5]."\t".$e[6]."\t".$e[7]."\t".$e[8]."\t".$e[9]."\t".$e[10];
	}
	else {
		print $_;
	}
}

close IN;

exit(0);
