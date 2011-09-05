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
#  Module: Parser::ML::GraphMapping::graphmap_remove_redundancy.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 ecoprior map.list\n";
my $eco_priority = shift or die $usage;
my $mapfile = shift or die $usage;

### Read ecotype priority list
my @eco_list = ();
open ECO, $eco_priority or die "Cannot open ecotype priority file\n";
while( <ECO> ) {
	chomp;
	push(@eco_list, $_);
}

open OUTFILE, ">$mapfile.nonredundant" or die "Cannot open outfile\n";
open FILE, $mapfile or die "Cannot open map file\n";

my %hit_correction = ();
my %last_lines = ();
my $compare_chr = -1;
my $compare_pos = -1;
my $compare_id  = -1;
my $compare_ori = "";

### First round: find redundant read IDs mapped at the same position of the genome due to multi ecotype hits
while (my $line = <FILE>) {
	my @e = split(" ", $line);

	### Case: entry redundant with last entry
	if( $e[1] == $compare_chr && $e[2] == $compare_pos && $e[6] == $compare_id && $e[7] eq $compare_ori ) {
		$last_lines{$e[0]} = $line;
	}

	### Case: entry not redundant with last entry
	else {

		### Print stored lines non redundant according to ecotype priority
		if( $compare_chr != -1 ) {

			# Store decrease of hits counter
			my $hit_decrease = scalar(keys %last_lines) - 1;
			if( $hit_decrease > 0 ) {
				$hit_correction{$compare_id} += $hit_decrease;
			}

			# Print entry for ecotype with highest priority
			for(my $i =0; $i < @eco_list; $i++) {
				if( exists $last_lines{$eco_list[$i]} ) {
					print OUTFILE $last_lines{$eco_list[$i]};
					last;
				}
			}
		}
		
		### Reset variables
		%last_lines = ();
		$last_lines{$e[0]} = $line;
		$compare_chr = $e[1];
		$compare_pos = $e[2];
		$compare_id  = $e[6];
		$compare_ori = $e[7];
	}
}
close FILE;
close OUTFILE;

### Second round: adapt hits counter
open OUTFILE, ">$mapfile.nonredundant.hitcorrected" or die "Cannot open outfile\n";
open FILE, "$mapfile.nonredundant" or die "Cannot open nonredundant map file\n";
while (my $line = <FILE>) {
	my @e = split(" ", $line);

	# Case: Hits have to be adapted
	if( exists $hit_correction{ $e[6] } ) {
		my $hits = $e[9] - $hit_correction{ $e[6] };

		print OUTFILE	$e[0] ."\t". $e[1] ."\t". $e[2] ."\t". $e[3] ."\t".$e[4] ."\t". $e[5] ."\t". $e[6] ."\t". $e[7] ."\t".
				$e[8] ."\t". $hits ."\t". $e[10] ."\t". $e[11] ."\t". $e[12] ."\t". $e[13] ."\t". $e[14] ."\t". $e[15] ."\n";

	}

	# Case: Hits stay the same
	else {
		print OUTFILE $line;
	}
}
close FILE;
close OUTFILE;

exit(0);
