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
#  Module: Parser::ML::select_chromosomes.pl
#  Purpose:
#  In:
#  Out:
#

my $usage    = "$0 chr_list stop file\n\n";
my $chr_list = shift or die "Please specify chromosome list\n\n$usage\n\n";
my $file     = shift or die "Please specify map.list file\n\n$usage\n\n";
my $stop     = shift;

if(! defined $stop) {$stop = -1};

my @chr_tmp = split(",", $chr_list);
my %chr = ();
foreach(@chr_tmp) { $chr{$_} = 1; }

open FILE, $file or die "cannnot open $file\n";

while(<FILE>) {
	my @entries = split(/\t/, $_);
	if( exists $chr{$entries[0]} ) {
		print $_;
	}
	if($entries[0] == $stop) { last; }
}

exit(0);
