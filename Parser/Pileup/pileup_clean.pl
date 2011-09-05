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
#  Module: Parser::Pileup::pileup_clean.pl
#  Purpose:
#  In:
#  Out:
#

###############################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger 
#Date 		07/04/07
#Version	0.1
#Input		Check pileup variant file for correctness
###############################################################


my $Usage = "\n\n$0 quality minlen minsup minfreq pileup_file\n\n";

my $quality = shift or die $Usage;
my $minlen  = shift or die $Usage;
my $minsup  = shift or die $Usage;
my $minfreq = shift or die $Usage;
my $file    = shift or die $Usage;

open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	if( $#a >= 10 ) {
		if( ($a[5] >= $quality) && (length($a[8]) >= $minlen) && ($a[10] >= $minsup) && ($a[10]/$a[7] >= $minfreq) ) {
			print $_;
		}
	}
}

close FILE;

exit(0);
