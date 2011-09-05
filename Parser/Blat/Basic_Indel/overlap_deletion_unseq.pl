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
#  Module: Parser::Blat::Basic_Indel::overlap_deletion_unseq.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 deletions unsequenced\n";


my $del = shift or die $usage;
open UNSEQ, shift or die $usage;

my %U = ();

while (<UNSEQ>) {
	my @a = split " ";
	for (my $i = $a[2]; $i < $a[3]; $i++) {
		$U{$a[1]."#".$i} = 1;
	}
}

close UNSEQ;

open FILE, $del or die "Cannot open file\n";

while (<FILE>) {
	chomp;
	my @a = split " ";
	my $count = 0;
	
	for (my $i = $a[2]; $i <= $a[3]; $i++) {
		if (defined($U{$a[1]."#".$i})) {
	                $count++;
		}
        }

	print $_, "\t", $count, "\n";
}

