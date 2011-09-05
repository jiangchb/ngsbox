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
#  Module: Parser::FL::parse_reads_by_id.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 idfile readsfile\n";
my $idfile = shift or die $usage;
my $readfile = shift or die $usage;

open IDFILE, $idfile or die $usage;
open READFILE, $readfile or die $usage;

my %IDS = ();

while (my $line = <IDFILE>) {
	chomp($line);
	$IDS{$line} = 1;
}

close IDFILE;

while (my $line = <READFILE>) {
	my @a = split " ", $line;
	if (defined($IDS{$a[0]})) {
		print $line;
	}
}

close READFILE;

