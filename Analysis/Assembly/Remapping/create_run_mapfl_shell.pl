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
#  Module: Analysis::Assembly::Remapping::create_run_mapfl_shell.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 project_folder\n";

my $project_folder = shift or die $usage;




my @strain_folders = glob($project_folder."/*");

# Parse strain folders
my @strain_folders = glob($project_folder."/*");
foreach my $strain_folder (sort @strain_folders) {

	# Parse run folders
	my @run_folders = glob($strain_folder."/run_*");
	foreach my $run_folder (sort @run_folders) {

		# ~/shore_mapfl/shore mapflowcell -i /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked_All_Mapper/TAIR8.v1.Masked.fa.shore -f . -u -n 7% -g 1 -D 50000 -S -L 14 -R 1 -c 16 -b 200000 -v blat -p
		print "~/shore_mapfl/shore mapflowcell -i /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked_All_Mapper/TAIR8.v1.Masked.fa.shore -f $run_folder -u -n 7% -g 1 -D 1000 -S -L 14 -R 1 -c 16 -b 200000 -v blat -p\n";
	}
}

