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
#  Module: Analysis::SV::compare.pm
#  Purpose:
#  In:
#  Out:
#


package compare;

sub new {
        my $self = {
        };
        bless $self;
        return $self;
}

sub comp {
	my ($self, $min, $max, $chr1, $pos1, $end1, $chr2, $pos2, $end2) = @_;

	my $id_var = ($max - $min + 1) / 2;

	## Chromosome need to be the same
	if ($chr1 != $chr2) {
		return 0;
	}

	## Length of deletions need to similar
	my $len1 = $end1-$pos1+1;
	my $len2 = $end2-$pos2+1;
	my $len_diff = abs($len1-$len2);

	if ($len_diff > $id_var) {
		return 0;
	}

	## Overlaps need to large enough to leave space for the deletion
	my $comb_pos = $pos1 > $pos2 ? $pos1 : $pos2;
	my $comb_end = $end1 < $end2 ? $end1 : $end2;	
	my $comb_len = $comb_end - $comb_pos + 1;

	if ($comb_len < 1 or $comb_len < ($len1 - $id_var) or $comb_len < ($len2 - $id_var)) {
		return 0;
	}

	return 1;

}



1;

