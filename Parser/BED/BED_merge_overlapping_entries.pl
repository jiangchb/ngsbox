#!/usr/bin/perl
use strict;
use warnings;

my $extension = shift;# or die;
my $file = shift;# or die;
open IN, $file or die "Cannot open input file\n";

my $last_chr = "NA";
my $last_start = -999999;
my $last_end = -999999;


while( <IN> ) {

	chomp;
	my ($chr, $beg, $end, $name) = split(/\t/, $_);

	# Extend region
	$beg -= $extension;
	$end += $extension;
	if($beg < 1) { $beg = 1; }

	if($chr ne $last_chr) {
		if($last_chr ne "NA") {
			print "$last_chr\t$last_start\t$last_end\tenriched\n";
		}
		if($beg > 1) {
			my $depletion_end = $beg - 1;
			print "$chr\t1\t$depletion_end\tdepleted\n";
		}
		
		$last_chr = $chr;
		$last_start = $beg;
		$last_end = $end;
	}
	elsif($beg <= $last_end) {
		$last_end = $end;
	}
	else {
		my $depletion_start = $last_end + 1;
		my $depletion_end = $beg - 1;
		print "$chr\t$last_start\t$last_end\tenriched\n";
		print "$chr\t$depletion_start\t$depletion_end\tdepleted\n";

		$last_start = $beg;
		$last_end = $end;
	}
}

close IN;

exit(0);
