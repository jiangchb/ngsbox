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
#  Module: Parser::SOAP::parse_gapped_reads_gm.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 gmoutput\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;


my %IDS = ();
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if ($a[2] =~ m/\-/) {
		if (not defined($IDS{$a[3]})) {
			print $a[3], "\n";
			$IDS{$a[3]} = 1;
		}
	}
}


