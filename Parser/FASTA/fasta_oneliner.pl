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
#  Module: Parser::FASTA::fasta_oneliner.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 fastafile\n";

my $in = shift or die $usage;

open IN, $in or die $usage;

my $seq = "";
my $id  = "";

while (my $line = <IN>) {
	chomp($line);

	if (substr($line, 0, 1) eq ">") {

		if ($seq ne "") {
			print "$id\n$seq\n";
		}
		$seq = "";
		$id = $line;
	}
	else {
		$seq .= $line;
	}
}

if ($seq ne "") {
	print "$id\n$seq\n";
}

	

