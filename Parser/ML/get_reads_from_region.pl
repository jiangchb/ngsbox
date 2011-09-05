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
#  Module: Parser::ML::get_reads_from_region.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 map.list regions readlength \n";
my $map_list = shift or die $usage;
my $regions = shift or die $usage;
my $read_length = shift or die $usage;

if(! defined $read_length) { $read_length = 0; }

### Read target regions
my %REGIONS = ();
open REG, $regions or die "cannot open $regions\n";
while(<REG>) {
	my @a = split " ", $_;
	for (my $i = $a[1]-$read_length; $i<=$a[2]; $i++) {
		$REGIONS{$a[0]."#".$i} = 1;
	}
}
close REG;

### Get reads from target region
open FILE, $map_list or die "cannnot open $map_list\n";
while (<FILE>) {
	my @a = split " ", $_;
	if (defined($REGIONS{$a[0]."#".$a[1]})) {
		print $_;

		### Clean alignment from brackets and read-base
		my $seq = "";
		for (my $i = 0; $i<length($a[2]); $i++) {
			if (substr($a[2], $i, 1) eq "[") {
				$seq .= substr($a[2], $i+2, 1);
				$i+=3;
			}
			else {
				$seq .= substr($a[2], $i, 1);
			}
		}
		print ">$a[3]\n$seq\n";
	}
}


exit(0);
