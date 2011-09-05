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
#  Module: Parser::ML::select_for_assembly.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 maxhit file\n";

my $max_hit = shift or die $usage;
my $file    = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while ( <IN> ) {
	my @elements = split(/\t/, $_);

	if(@elements != 13) {
		print STDERR "Not enough columns: $elements[0]\t$elements[1]\n";
	}
	elsif( $elements[6] !~ /\d+/ ) {
		print STDERR "Hits not int: $elements[0]\t$elements[1]\t$elements[6]\n";
	}
	elsif(	($elements[6] == 1) || ($elements[6] <= $max_hit && $elements[5] == 0) ) {
		print $_;
	}
}

close IN;

exit(0);
