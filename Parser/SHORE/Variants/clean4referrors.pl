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
#  Module: Parser::SHORE::Variants::clean4referrors.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 quality_variants.txt referrors.txt\n";

my $qual_file = shift or die $usage;
my $err_file = shift or die $usage;

my %ERR = ();

open FILE, $err_file or die "Cannot open file\n";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	$ERR{$a[0]."#".$a[1]} = 1;
}

close FILE;

open FILE, $qual_file or die "Cannot open file\n";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (not defined($ERR{$a[1]."#".$a[2]})) {
		print $line;
	}
}


