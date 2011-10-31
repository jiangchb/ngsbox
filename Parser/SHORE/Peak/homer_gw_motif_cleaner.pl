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
#  Module: Parser::SHORE::Peak::
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 infile\n\n";
my $file = shift or die $usage;

my %hash = ();


### Print distinct binding sites (distinct on chr, pos)
open FILE, $file or die $usage;
while (<FILE>) {

	chomp;
	my @a = split(/\t/, $_);

	my $chr_pos = $a[1] . "#" . $a[2];

	if(exists $hash{$chr_pos} ) {

		my @b = split(/\t/, $hash{$chr_pos});

		if($a[5] > $b[5]) {

			$hash{$chr_pos} = $_;
		}
	}
	else {
		$hash{$chr_pos} = $_;
	}
}
close FILE;


foreach my $entry (sort keys %hash) {
	print $hash{$entry} . "\n";
}

exit(0);

