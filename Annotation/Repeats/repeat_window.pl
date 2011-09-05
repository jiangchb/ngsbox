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
#  Module: Annotation::Repeats::repeat_window.pl
#  Purpose:
#  In:
#  Out:
#



my $start = shift;
my $win   = shift;
my $file = shift;

open FILE, $file or die;
my $last_len = $start;


### Plot count in sliding window
#my $win_count = 0;
#while( <FILE> ) {
#	my ($len, $count, $sum) = split("\t", $_);
#
#	if($len >= $start) {
#		if( ($len%$win==0) && ($len > $start) ) {
#			print "$last_len\t$win_count\n";
#			$win_count = 0;
#			$last_len = $len;
#		}
#
#		$win_count+=$count;
#	}
#}


### Plot sum up to sliding window
my $new_sum = 0;
while( <FILE> ) {
	my ($len, $count, $sum) = split("\t", $_);

	if($len >= $start) {
		if( $len%$win == 0) {
			print "$len\t$new_sum\n";
		}
		$new_sum+=$count;
	}
}
