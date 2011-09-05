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
#  Module: Parser::FASTA::fasta_2_flat.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 <fastafile> <outfile prefix>\n";

my $in  = shift or die $usage;
my $out = shift or die $usage;

open IN, $in or die $usage;
open OUT1, ">$out.1" or die $usage;
open OUT2, ">$out.2" or die $usage;


while (my $header1 = <IN>) {

	my $seq1    = <IN>;
	my $header2 = <IN>;
	my $seq2    = <IN>;

	chomp($header1);
	chomp($header2);
	chomp($seq1);
	chomp($seq2);


	# Quality string
	my $q2 = "";

	for (my $i = 0; $i < length($seq1); $i++) {
		$q2 .= "I";
	}

	print OUT1 substr($header1, 1) . "\t$seq1\t1\t$q2\n";
	print OUT2 substr($header1, 1) . "\t$seq2\t2\t$q2\n";
}
