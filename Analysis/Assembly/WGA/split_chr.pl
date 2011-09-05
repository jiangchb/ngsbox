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
#  Module: Analysis::Assembly::WGA::split_chr.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 fastafile\n";

my $file = shift or die $usage;
open FILE, $file or die $usage;

my $seq  = "XX";
my $id   = "";
my $desc = "";
my $chr  = "";

### Parse assembly file
while (my $line = <FILE>) {
	chomp($line);

	# Header found
	if (substr($line, 0, 1) eq ">") {
		
		# Split header
		my @e = split(" ", $line);
		my @f = split("_", $e[2]); 

		if ($seq ne "XX") {
			print OUT "$id | $desc\n$seq\n";

			# Open new file if new chromosome arm is found
			if($f[0] ne $chr) {
				close OUT;
				$chr = $f[0];
				open OUT, ">$chr.fa" or die "Cannot open $chr.fa\n";
			}
		}
		else {
			$chr = $f[0];
			open OUT, ">$chr.fa" or die "Cannot open $chr.fa\n";
		}

		# Reset
		$seq  = "";
		$id   = $e[0];
		$desc = $e[2];
		$chr  = $f[0];
	}
	else {
		$seq .= $line;
	}
}

### Last one
if ($seq ne "") {
	print OUT "$id | $desc\n$seq\n";
}

	

