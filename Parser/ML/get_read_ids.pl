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
#  Module: Parser::ML::get_read_ids.pl
#  Purpose:
#  In:
#  Out:
#


my $map_list = shift;
my $id_list = shift;

### Store IDs
my %IDS = ();
open IDLIST, $id_list or die "cannot open $id_list\n";
while(<IDLIST>) {
	chomp;
	$IDS{$_} = 1;
}
close IDLIST;

### Get alignments with specified IDs
open FILE, $map_list or die "cannnot open $map_list\n";
my $count = 0;
while (<FILE>) {
	my @a = split " ", $_;

	### DEBUG #######################
	if($count > 100000) {
		print STDERR "$a[0]\t$a[1]\n";
		$count = 0;
	}
	$count++;
	#################################
	
	if ( exists $IDS{$a[3]} ) {
		print $_;
	}
}

exit(0);
