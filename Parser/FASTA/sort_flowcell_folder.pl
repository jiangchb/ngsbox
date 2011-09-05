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
#  Module: Parser::FASTA::sort_flowcell_folder.pl
#  Purpose:
#  In:
#  Out:
#




### User params
my $usage = "$0 lane\n";

my $lane         = shift or die $usage;

### Get all subfolders (AMOS batches)
my @subfolders = glob("$lane/*/*");


### Debug output
foreach my $subfolder (@subfolders) {
	print "$lane/$subfolder/reads_0.fl\n";
}

foreach my $subfolder (@subfolders) {
	system("sort -n -k1 $subfolder/reads_0.fl > $subfolder/reads_0.tmp");
	system("mv $subfolder/reads_0.tmp $subfolder/reads_0.fl");
}

exit(0);
