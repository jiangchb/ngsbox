#!/usr/bin/perl
use strict;
use warnings;

my $extension = shift;# or die;
my $file = shift or die;
open IN, $file or die "Cannot open input file\n";

my $last_chr = "NA";
my $last_beg = -999999;
my $last_end = -999999;


while( <IN> ) {

	chomp;
	my ($chr, $beg, $end, $name) = split(/\t/, $_);

	if($chr ne $last_chr || $beg != $last_beg || $end != $last_end) { 
		$beg -= $extension;
		$end += $extension;
		if($beg < 1) { $beg = 1; }

		print "$chr\t$beg\t$end\t$name\n";
	}
	
	$last_chr = $chr;
	$last_beg = $beg;
	$last_end = $end;
}

close IN;

exit(0);
