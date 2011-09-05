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
#  Module: Parser::FASTA::fasta_2_batch.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 file batchnumber\n";
my $file = shift or die $usage;
my $batch = shift or die $usage;
if ($batch < 1) {
	die ("batchnumber must be larger than 1\n");
}

my $num = `grep -c '^>' $file`;

my $size = ($num / $batch) + 1;

open FILE, $file or die "Cannot open $file\n";

my $count = 0;
my $filenum = 1;
open OUT, ">".$file.".".$filenum or die "Cannot open out file\n";
while (<FILE>) {
	if (substr($_, 0, 1)  eq ">") {
		$count++;
		if ($count >= $size) {
			$filenum++;
			close OUT or die "Cannot close out file\n";
			open OUT, ">".$file.".".$filenum or die "Cannot open out file\n";
			$count = 1;
		}
	}
	print OUT $_;
}



exit(0);
