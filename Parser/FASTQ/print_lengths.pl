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
#  Module: Parser::FASTA::print_lengths.pl
#  Purpose:
#  In:
#  Out:
#

my $usage="$0 fasta\n";
open FILE, shift or die $usage;

my $seq = "";

while (<FILE>) {
	if (substr($_, 0, 1) eq ">") {
		print length($seq), "\n" if ($seq ne "");
		$seq = "";
	}
	else {
		chomp();
		$seq.=$_;
	}
}
print length($seq), "\n" if ($seq ne "");



