#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Parser::Pileup::evaluate_strand_bias.pl
#  Purpose:
#  In:
#  Out:
#



my $pileup = shift or die "Please specify pileup file\n";

open PILEUP, $pileup or die "Cannot open input file\n";

# Loop through pileup
while( <PILEUP> ) {
	chomp;

	my @a = split("\t", $_);

	if( ($a[7] <= 1) && ($a[6] >= 10) ) {

		my $fwd = 0;
		my $rev = 0;

		for(my $i = 0; $i < length($a[9]); $i++) {

			if(substr($a[9], $i, 1) eq ".") {
				$fwd++;
			}
			elsif(substr($a[9], $i, 1) eq ",") {
				$rev++;
			}
		}

		#if( ($fwd + $rev) >= 10) {
		#	my $ratio = &min($fwd, $rev) / &max($fwd, $rev);
		#	print "$a[0]\t$a[1]\t$ratio\n";
		#}

		if( ($fwd + $rev) >= 10) {
			my $ratio = 0;
			
			if($rev > $fwd) {
				$ratio = 2 - ($fwd / $rev);
			}
			else {
				$ratio = $rev / $fwd;
			}

			print "$a[0]\t$a[1]\t$ratio\n";
		}
	}
}

close PILEUP;

exit(0);

sub min {
	my ($a, $b) = @_;
	return $a if $a < $b;
	return $b;
}
sub max {
	my ($a, $b) = @_;
	return $a if $a > $b;
	return $b;
}

