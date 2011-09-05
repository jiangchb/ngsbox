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
#  Module: Analysis::Filter::low_complexity_filter.pl
#  Purpose:
#  In:
#  Out:
#



my $file = shift;

open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;
	my @elements = split(/\t/, $_);
	if( &seq_complexity($elements[1]) ) {
		print "$_\n";
	}
}
close IN;

exit(0);


sub seq_complexity
{
        my $seq = shift;
	my $length = length($seq);
        my %occ   = (A => 0, C => 0,G => 0,T => 0);
	my %count = (A => 0, C => 0,G => 0,T => 0);

        for(my $i = 0; $i < length($seq); $i++) {
                $count{substr($seq, $i, 1)}++;
		$occ{substr($seq, $i, 1)} = 1;
        }

        if( 
		( ($occ{A} + $occ{C} + $occ{G} + $occ{T}) >= 3) &&
		( ($length - $count{A}) > 4 ) &&
		( ($length - $count{C}) > 4 ) &&
		( ($length - $count{G}) > 4 ) &&
		( ($length - $count{T}) > 4 )
	) {
                return(1);
        }
        else {
                return(0);
        }
}
