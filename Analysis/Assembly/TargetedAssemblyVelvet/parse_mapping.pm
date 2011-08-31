#!/usr/bin/perl
###############################################################################
# Author 	Korbinian Schneeberger, Stephan Ossowski 
# Date 		09/26/07
# Version	0.3
# Input		map.list, chr, start, end
# Function	returns all reads starting in a specified region of the mapping
###############################################################################

use strict;
use warnings;

package parse_mapping;
 
sub get
{
	my ($file, $chr, $begin, $end) = @_;
	my @results = ();
	open FILE, $file or die "Cannot open $file\n";
	open IDX, "$file.idx" or die "Cannot open $file.idx\n";

	### Read in index file an jump FILE to right position
	my $jump = 0;
	SETTING: while ( <IDX> ) {
		my @a = split(" ", $_);
		if ($a[0] > $chr or ($a[0] == $chr and $a[1] >= $begin)) {
			last SETTING;
		}
		else { $jump = $a[2]; }
	}
	seek(FILE, $jump, 0);
	
	### Set file to region of interest and read in first fragment
	FIND: while (<FILE>) {
		my ($current_chr, $current_pos) = split (" ", $_);
		if ($current_chr == $chr and $current_pos>= $begin) {
			chomp;
			push @results, $_;
			last FIND;
		}
	}

	### Parse until the end of region of interest
	REG: while (<FILE>) {
		chomp;
		my ($current_chr, $current_pos) = split (" ", $_);
		if ($chr != $current_chr or $current_pos > $end) {
			last REG;
		}
		push @results, $_;
	}

	return(@results);
}
1;
