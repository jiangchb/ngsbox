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
#  Module: Support::qsub::qsub_import_fastq.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 infolder outfolder qsubName\n\n";
my $infolder  = shift or die $usage;
my $outfolder = shift or die $usage;
my $qsub_name = shift or die $usage;

open OUT, ">$qsub_name" or die $usage;


### Print qsub header and exports
print OUT "#!/bin/bash\n\n";
print OUT "#\$ -e $outfolder\n";
print OUT "#\$ -o $outfolder\n\n";
print OUT "export SHORE=/users/GD/tools/shore/\n";
print OUT "export TMPDIR=/users/GD/projects/familiaMarruecos/tmp/\n\n";

my $counter = 1001;

my @subfolders = glob($infolder . "/*");

foreach my $subfolder (@subfolders) {

	my @subpath = split("/", $subfolder);
	my $subleaf = $subpath[$#subpath];

	my @files = glob($subfolder . "/*");

	my @read1 = ();
	my @read2 = ();

	foreach my $file (@files) {
		my @filepath = split("/", $file);
		my $fileleaf = $filepath[$#filepath];

		my @junk = split("_", $fileleaf);
	
		if($junk[2] == 1) {
			push @read1, $file;
		}
		else {
			push @read2, $file;
		}
	}

	print OUT "/users/GD/tools/shore/shore import -v Fastq -e Shore -a genomic -i $counter -o $outfolder/$subleaf -n 4 -g -c -V 10 -x ";
	print OUT join ",", @read1;
	print OUT " -y ";
	print OUT join ",", @read2;
	print OUT "\n\n";

	$counter++;
}
