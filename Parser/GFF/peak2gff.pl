#!/usr/bin/perl

###########################################################
# Create gff for coverage by position plot
# Author: Stephan Ossowski
# Date: 05/03/2008
###########################################################

use strict;
use warnings;

my $assay         = shift;
my $coverage_file = shift;
my $enriched_file = shift;

### Parse Segments
my $zero_start = 1;
my %seg_hash = ();

open EF, $enriched_file or die "Cannot open segment file\n";
while( <EF> ) {
	my @elem = split("\t", $_);
	my $zero_end = $elem[1] - 1;

	# print entry for non-enriched region
	print "$elem[0]\tIlluminaGA2\t$assay\t$zero_start\t$zero_end\t0\t.\t.\tName=$assay\n";

	# Store enriched segment
	for(my $x = $elem[1]; $x < $elem[2]; $x++) {
		$seg_hash{$elem[0] ."-". $x} = 1;
	}
	$zero_start = $elem[2] + 1;
}
close(EF);


### Parse coverage
my $chr   = -1;
my $start = -1;
my $end   = -1;
my $cov   = -1;

open COV, $coverage_file or die "Cannot open coverage file\n";

while( <COV> ) {
	my @elem = split("\t", $_);

	if( ($elem[2] != $cov) || ($elem[0] ne $chr) || ($elem[1] > $end + 1) ) {

		# Print enriched region
		if( (exists $seg_hash{ $chr . "-" . $start }) && ($chr != -1) ) {
			print "$chr\tIlluminaGA2\t$assay\t$start\t$end\t$cov\t.\t.\tName=$assay\n";
		}

		# reset
		$chr = $elem[0];
		$start = $elem[1];
		$end = $elem[1];
		$cov = $elem[2];
	}

	$end = $elem[1];
}
close(COV);

exit(0);

