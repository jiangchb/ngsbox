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
#  Module: Parser::ML::fix_position.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 file\n\n";
my $file  = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	my @e = split(/\t/, $_);

	if( $e[2] =~ /\[L/ ) {
		$e[1]++;
		print $e[0]."\t".$e[1]."\t".$e[2]."\t".$e[3]."\t".$e[4]."\t".$e[5]."\t".$e[6]."\t".$e[7]."\t".$e[8]."\t".$e[9]."\t".$e[10];
	}
	else {
		print $_;
	}
}

close IN;

exit(0);
