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
#  Module: Parser::FASTQ::split_fastq.pl
#  Purpose:
#  In:
#  Out:
#



my $split = shift;
my $file  = shift;

my $count = 0;
my $file_count = 1;

open IN, $file or die "Cannot open input file\n";
open OUT, ">$file.$file_count" or die "Cannot open output file\n";


my $line = "";
while( $line = <IN> ) {
	# Print 4 lines of fastq entry
	print OUT $line;
	$line = <IN>;
	print OUT $line;
	$line = <IN>;
	print OUT $line;
	$line = <IN>;
	print OUT $line;

	$count++;

	if($count >= $split) {
		if($file_count > 2) { exit; }

		$count = 0;
		$file_count++;
		close OUT;
		open OUT, ">$file.$file_count" or die "Cannot open output file\n";
	}
}
close IN; close OUT;

exit(0);
