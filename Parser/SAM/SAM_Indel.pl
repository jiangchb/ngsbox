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
#  Module: Parser::SAM::SAM_Indel.pl
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


my $file = shift or die "\n\nUsage: $0 SAMfile\n\n";
open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	#my $pe = 0;
	#if ( hex($a[1]) & 0x0040 ) {
	#	$pe = 1;
	#}
	#elsif ( hex($a[1]) & 0x0080 ) {
	#	$pe = 2;
	#}
	
	if( ($a[5] =~ /D/) || ($a[5] =~ /I/) ) {
		my @c = split(/[MIDS]/, $a[5]);
		
		if( $c[0] >= 8 && $c[$#c] >= 8) {
			print $_;
		}
	}
}

close FILE;

exit(0);
