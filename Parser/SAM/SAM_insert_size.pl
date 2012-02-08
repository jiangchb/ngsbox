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


my $usage = "\n\nUsage: $0 min max SAMfile\n\n";

my $min = shift or die $usage;
my $max = shift or die $usage;
my $file = shift or die $usage;


open FILE, $file;

while(<FILE>) {
	if(substr($_, 0, 1) eq "@") {
		print $_;
	}
	else {
		my @a = split("\t", $_);

		my $in_range = 0;

		if( abs($a[8]) >= $min && abs($a[8]) <= $max) {
			$in_range = 1;
		}

		if($in_range == 1) {
			print $_;
		}
	}
}

close FILE;

exit(0);
