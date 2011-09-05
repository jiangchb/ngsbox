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
#  Module: Parser::FASTQ::merge_read_1_2.pl
#  Purpose:
#  In:
#  Out:
#



my $usage  = "$0 file1 file2\n";
my $file1   = shift or die $usage;
my $file2   = shift or die $usage;

open F1, $file1 or die "Cannot open input read 1 file\n";
open F2, $file2 or die "Cannot open input read 2 file\n";

while( <F1> ) {
	my $sh1   = $_;
	my $seq1  = <F1>;
	my $qh1   = <F1>;
	my $qual1 = <F1>;

	my $sh2   = <F2>;
	my $seq2  = <F2>;
	my $qh2   = <F2>;
	my $qual2 = <F2>;

	print "$sh1$seq1$qh1$qual1$sh2$seq2$qh2$qual2";
}

close F1;
close F2;

exit(0);
