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
#  Module: Support::qsub::qsub_merge_maplist.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 runfolder outfolder qsubname\n\n";
my $runfolder  = shift or die $usage;
my $qsub_name  = shift or die $usage;

open OUT, ">$qsub_name" or die $usage;


### Print qsub header and exports
print OUT "#!/bin/bash\n\n";
print OUT "#\$ -e $runfolder\n";
print OUT "#\$ -o $runfolder\n\n";
print OUT "source /users/GD/so/sossowski/.bashrc\n";
print OUT "export TMPDIR=/users/GD/projects/familiaMarruecos/tmp/\n\n";


my $counter = 1001;

my @lanefolders = glob($runfolder . "/*");

foreach my $lanefolder (@lanefolders) {

	my @typefolders = glob($lanefolder . "/*");

	foreach my $typefolder (@typefolders) {
		
		my @subpath = split("/", $typefolder);
		my $subleaf = $subpath[$#subpath];

		my @lengthfolders = glob($typefolder . "/*");

		foreach my $lengthfolder (@lengthfolders) {
			print OUT "mv $lengthfolder/map.list.gz $lengthfolder/map.first.gz\n\n";
			print OUT "mv $lengthfolder/map.list.gzx.gz $lengthfolder/map.first.gzx.gz\n\n";

			if($subleaf ne "single") {
				print OUT "shore sort -m -k4i -o $lengthfolder/map.list.gz $lengthfolder/map.first.gz $lengthfolder/map.left.gz\n\n";
			}
			else {
				print OUT "shore sort -m -k1i2i -o $lengthfolder/map.list.gz $lengthfolder/map.first.gz $lengthfolder/map.left.gz\n\n\n\n";
			}
		}
	}
}
