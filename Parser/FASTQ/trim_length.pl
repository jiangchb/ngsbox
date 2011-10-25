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
#  Module: Parser::FASTQ::cut_length.pl
#  Purpose:
#  In:
#  Out:
#



my $usage  = "$0 start end file\n";
my $beg    = shift or die $usage;
my $end    = shift or die $usage;
my $file   = shift or die $usage;

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	my $h1  = $_;
	my $seq = <IN>;
	my $h2  = <IN>;
	my $qual = <IN>;

	chomp $seq;
	chomp $qual;
	my $print_seq  = substr( $seq,  $beg - 1, ($end - $beg + 1) );
	my $print_qual = substr( $qual, $beg - 1, ($end - $beg + 1) );

	print "$h1$print_seq\n$h2$print_qual\n";
}

close IN;

exit(0);
