#! /usr/bin/perl
use strict;
use warnings;

my $usage = "\n$0 ngsbox_folder\n\n";
my $infolder  = shift or die $usage;



### Dir Level 1
my @l1_folders = glob($infolder . "/*");

foreach my $l1_folder (@l1_folders) {


	### Dir Level 2
	my @l1_path = split("/", $l1_folder);
	my $l1_leaf = $l1_path[$#l1_path];

	my @l2_folders = glob($l1_folder . "/*");


	foreach my $l2_folder (@l2_folders) {
		

		### Dir Level 3
		my @l2_path = split("/", $l2_folder);
		my $l2_leaf = $l2_path[$#l2_path];


		my @l3_files = glob($l2_folder . "/*");

		foreach my $l3_file (@l3_files) {

			my @l3_path = split("/", $l3_file);
			my $l3_leaf = $l3_path[$#l3_path];
	
			### Only modify text files
			if( (-T $l3_file) && ($l3_leaf =~ /\.pl/) ) {


my @mod_file = ("#! /usr/bin/perl
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

\n");

				open IN, $l3_file or die "\nError: Cannot open original perl script\n\n";
				while( <IN> ) {
					if( ($_ !~ /\/usr\/bin\/perl/) && ($_ !~ /use strict;/) && ($_ !~ /use warnings;/) ) {
						push(@mod_file, $_);
					}
				}
				close IN or die;

				open OUT, ">$l3_file" or die "\nError: Cannot open modified perl script\n\n";

				print OUT @mod_file;
			}
		}
	}
}
