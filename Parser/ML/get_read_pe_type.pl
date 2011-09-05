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
#  Module: Parser::ML::get_read_pe_type.pl
#  Purpose:
#  In:
#  Out:
#


my $map_list = shift;
my $pe_list = shift;

### Get PE types
my %pe_types = ();
my @tmp = split(",", $pe_list);
foreach ( @tmp ) { $pe_types{$_} = 1; }

### Get alignments with specified PE types
open FILE, $map_list or die "cannnot open $map_list\n";
while (<FILE>) {
	my @a = split("\t", $_);
	if ( exists $pe_types{$a[9]} ) {
		print $_;
	}
}

exit(0);
