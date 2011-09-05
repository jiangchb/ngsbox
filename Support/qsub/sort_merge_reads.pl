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
#  Module: Support::qsub::sort_merge_reads.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 infolders\n\n";
my $ins  = shift or die $usage;


my @infolders = split(",", $ins);
my $files_r1 = "";
my $files_r2 = "";


foreach my $infolder (@infolders) {

	### Dir Level 1
	my @typefolders = glob($infolder . "/*");

	foreach my $typefolder (@typefolders) {

		### Dir Level 2
		my @typepath = split("/", $typefolder);
		my $typeleaf = $typepath[$#typepath];

		my @lengthfolders = glob($typefolder . "/*");

		foreach my $lengthfolder (@lengthfolders) {
		
			### Dir Level 3
			my @lengthpath = split("/", $lengthfolder);
			my $lengthleaf = $lengthpath[$#lengthpath];

			if($typeleaf eq "1") {
				if(-e "$lengthfolder/reads_0.fl") {
					$files_r1 .= "$lengthfolder/reads_0.fl,"
				}
				elsif(-e "$lengthfolder/reads_0.fl.gz") {
					$files_r1 .= "$lengthfolder/reads_0.fl.gz,"
				}
			}
			elsif($typeleaf eq "2") {
				if(-e "$lengthfolder/reads_0.fl") {
					$files_r2 .= "$lengthfolder/reads_0.fl,"
				}
				elsif(-e "$lengthfolder/reads_0.fl.gz") {
					$files_r2 .= "$lengthfolder/reads_0.fl.gz,"
				}
			}
		}	
	}
}

chop($files_r1);
chop($files_r2);
print "\n\nshore sort -o LongIndels/reads_1.fl -m -k 1i -i $files_r1\n\n";
print "\n\nshore sort -o LongIndels/reads_2.fl -m -k 1i -i $files_r2\n\n";

