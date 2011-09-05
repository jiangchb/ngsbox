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
#  Module: Parser::Wiggle::coverage_distribution_in_wiggle.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 file\n";
my $file   = shift or die $usage;


### Read coverage idstribution from wiggle file
my %cov_dist = ();
open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;

	if( $_ =~ /track/ ) {
		# nothing for now
	}
	elsif( $_ =~ /fixedStep/ ) {
		my @a = split("\t", $_);
		# TODO: read chr, start, step, span
	}
	else {
		$cov_dist{$_}++;
	}
}
close IN;


### Print results
my $covsum = 0;
print "coverage\twindows\tpositions\tCumulative positions >= cov\n";
foreach my $cov (sort {$b<=>$a} keys %cov_dist) {

	my $positions = $cov_dist{$cov} * 20;
	$covsum += $positions;
	print "$cov\t" . $cov_dist{$cov} . "\t$positions\t$covsum\n";
}

exit(0);
