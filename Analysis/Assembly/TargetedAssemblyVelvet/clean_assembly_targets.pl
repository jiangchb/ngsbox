#!/usr/bin/perl
use strict;
use warnings;

my $flanking_region_length = shift;
my $region_file  = shift;

open IN, "$region_file" or die "Cannot open pc region file\n";

my $last_chr = -1;
my $last_beg = -10000;
my $last_end = -10000;

while(<IN>) {
	chomp;
	my ($chr, $beg, $end) = split (/\t/, $_);

	if( ( $beg > ($last_end + $flanking_region_length) ) || ( $chr != $last_chr) ) {

		my $region_beg = $last_beg - $flanking_region_length;
		my $region_end = $last_end + $flanking_region_length;
		if( $region_beg < 0 ) { $region_beg = 0; }

		if( ( ($region_end - $region_beg < 50000) ) && ($last_chr != -1) ) {
			print "$last_chr\t$region_beg\t$region_end\n";
		}
		$last_beg = $beg;
	}

	if( ($end >= $last_end) || ($chr != $last_chr) ) {
		$last_end = $end;
	}

	$last_chr = $chr;
}

my $region_beg = $last_beg - $flanking_region_length;
my $region_end = $last_end + $flanking_region_length;
if( $region_beg < 0 ) { $region_beg = 0; }

if($region_end - $region_beg < 50000) {
	print "$last_chr\t$region_beg\t$region_end\n";
}

exit(0);
