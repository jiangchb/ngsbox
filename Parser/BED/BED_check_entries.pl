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
#  Module: Parser::BED::BED_check_entries.pl
#  Purpose:
#  In:
#  Out:
#


my $file = shift;

my $last_chr = "";
my %cds = ();

open IN, $file or die "Cannot open input file\n";


while( <IN> ) {
	chomp;
	my ($chr, $beg, $end, $strand) = split(/\t/, $_);

	# New chromosome: reset CDS hash
	if($chr ne $last_chr) {
		$last_chr = $chr;
		%cds = ();
	}

	if($beg <= $end) {
		if( ! exists $cds{$beg} ) {
			print "$_\n";
			$cds{$beg} = $end;
		}
		elsif( $cds{$beg} != $end ) {
			print "$_\n";
			$cds{$beg} = $end;
		}
	}
	else {
		print STDERR "$_\n";
	}
}

close IN;

exit(0);
