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
#  Module: Parser::FASTA::fl_2_fasta.pl
#  Purpose:
#  In:
#  Out:
#


###############################################################
# Convert Shore read files into fasta format
# Authors 	Stephan Ossowski, Korbinian Schneeberger 
# Date 		07/04/07
# Version	?
# Input		Shore sequence read file (fl) 
# Output	FASTA file to STDOUT
###############################################################


my $file = shift;
open FILE, $file;

while(<FILE>) {
	my @a = split " ";
	print ">".$a[0]."\n".$a[1]."\n";

}

exit(0);
