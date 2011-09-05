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
#  Module: Parser::SAM::SAM_2_fl.pl
#  Purpose:
#  In:
#  Out:
#

###############################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger 
#Date 		07/04/07
#Version	0.1
#Input		Read flat sequence file and write fasta 
#		formated file to stdout
###############################################################


my $file = shift;
open FILE, $file;

open READ1, ">read1.fl" or die;
open READ2, ">read2.fl" or die;

while(<FILE>) {
	my @a = split("\t", $_);

	my $pe = 0;
	if ( hex($a[1]) & 0x0040 ) {
		$pe = 1;
	}
	elsif ( hex($a[1]) & 0x0080 ) {
		$pe = 2;
	}

	my $seq = "";
	my $qual = "";

	if( hex($a[1]) & 0x0010 ) {

		for(my $i = length($a[9]) - 1; $i >= 0; $i--) {

			# Reverse complement sequence
			$seq .= substr($a[9], $i, 1);

			# Reverse quality
			$qual .= substr($a[10], $i, 1);
		}

		$seq =~ tr/ACGTacgt/TGCAtgca/;
	}
	else {
		$seq = $a[9];
		$qual = $a[10];
	}

	if ( hex($a[1]) & 0x0040 ) {
		print READ1 $a[0] ."\t". $seq . "\t". $pe ."\t". $qual ."\n";
	}
	elsif ( hex($a[1]) & 0x0080 ) {
		print READ2 $a[0] ."\t". $seq . "\t". $pe ."\t". $qual ."\n";
	}
}

close FILE;
close READ1;
close READ2;

exit(0);
