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
#  Module: Parser::FL::parse_mapped_reads.pl
#  Purpose:
#  In:
#  Out:
#

######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Read prb values
######################################################################################


my $original_file = shift;
my $not_mapped_file = shift;

open (ALL, $original_file) or die "cannot open $original_file\n";
open (REST, $not_mapped_file) or die "cannot open $not_mapped_file\n";

my $rest_line = <REST>;
my ($rest_id, $rest_seq) = split("\t", $rest_line);

# Loop through original reads and get prb entry
while( my $all_line = <ALL> ) {
	chomp($all_line);
	my ($all_id, $all_seq, $all_prb, $all_qCal, $all_chas) = split(/\t/, $all_line);

	if($rest_id eq $all_id) {
		$rest_line = <REST>;
		if(! defined $rest_line) { $rest_line = "END\tNA"; }
		($rest_id, $rest_seq) = split("\t", $rest_line);
	}
	else {
		print "$all_id\t$all_seq\t$all_prb\t$all_qCal\t$all_chas\n";
	}
}

close ALL; close REST;

exit(0)
