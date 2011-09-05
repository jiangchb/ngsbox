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
#  Module: Analysis::Assembly::WGA::concatenate_deletions.pl
#  Purpose:
#  In:
#  Out:
#



my $file = shift;
open FILE, $file or die;

my $start_pos = 0;
my $last_pos = -1;
my $last_chr = -1;
my $reference = "";
my @e = ();
my @le = ();

while( <FILE> ) {
	chomp;
	@e = split("\t", $_);

	# Deletion entry
	if($e[4] eq "-") {

		# Deletion finished
		if( ($e[2] != $last_pos + 1) || ($e[1] != $last_chr) ) {

			my $length = $last_pos - $start_pos + 1;

			if($last_pos != -1) {
				print 	"$le[0]\t$last_chr\t$start_pos\t$last_pos\t$length\t$reference\t$le[5]\t$le[6]\t$le[7]\t$le[8]\n";
			}
			$start_pos = $e[2];
			$reference = "";

			for(my $i = 0; $i < scalar(@e); $i++) {
				$le[$i] = $e[$i];
				
			}
		}

		# Update positions and deleted sequence
		$last_chr = $e[1];
		$last_pos = $e[2];
		$reference .= $e[3];
	}
}

# Print last entry
my $length = $last_pos - $start_pos + 1;
print   "$le[0]\t$last_chr\t$start_pos\t$last_pos\t$length\t$reference\t$le[5]\t$le[6]\t$le[7]\t$le[8]\n";

exit(0);
