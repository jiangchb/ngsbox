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
#  Module: Analysis::Assembly::Remapping::create_run_corr4pe_shell.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 project_folder\n";

my $project_folder = shift or die $usage;



my $count = 0;
my @strain_folders = glob($project_folder."/*");

# Parse strain folders
my @strain_folders = glob($project_folder."/*");
foreach my $strain_folder (sort @strain_folders) {


	# Parse run folders
	my @run_folders = glob($strain_folder."/run_*");
	foreach my $run_folder (sort @run_folders) {
		
		my @lane_folders = glob($run_folder."/*");
		foreach my $lane_folder (sort @lane_folders) {

			my @a = split("/", $lane_folder);
			my $lane = $a[$#a];

			if( $lane =~ /[1-9]/ && $lane !~ /_/) {
				open OUT, ">>run_corr4pe_$count.sh" or die;
				print OUT "~/shore_mapfl/shore correct4pe -l $lane_folder -x 200 -e 1 -p -u\n";
				close OUT;
			}
		}
	}
	$count++;
}

