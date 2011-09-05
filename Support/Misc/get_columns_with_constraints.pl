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
#  Module: Support::Misc::get_columns_with_constraints.pl
#  Purpose:
#  In:
#  Out:
#


# select columns (comma seperated) and one constraint on one column (comma seperated "column,constrains") to filter files


my $usage = "$0 selected_columns constraint file\n";
my $columns    = shift or die $usage;
my $constraint = shift or die $usage;
my $file       = shift or die $usage;

my @colum = split(",", $columns);
my ($ccolumn, $cvalue) = split(",", $constraint);

open IN, $file or die "Cannot open input file\n";

while( <IN> ) {
	chomp;

	my @e = split("\t", $_);

	if( $e[$ccolumn - 1] > $cvalue ) {

		my $out = "";

		for(my $i = 0; $i <= $#colum; $i++) {

			my $value = $e[$colum[$i]-1];
		
			$out .= $value . "\t";
		}

		chop($out);

		print $out . "\n";
	}
}

close IN;

exit(0);
