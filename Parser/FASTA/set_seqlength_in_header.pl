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
#  Module: Parser::FASTA::set_seqlength_in_header.pl
#  Purpose:
#  In:
#  Out:
#



### User params
my $contig_file = shift;
my $min_length  = shift;

my $id = -1;
my %ctg_seq = ();

open CONTIG, $contig_file or die "Cannot open $contig_file\n";

while(<CONTIG>) {
	chomp($_);

	if (substr($_, 0, 1) eq ">") {
		$id = substr($_, 1);
	}
	else {
		$ctg_seq{$id} .= $_;
	}
}
close CONTIG;

### Print validated contigs
foreach my $id (sort keys %ctg_seq) {
	my $len = length($ctg_seq{$id});

	if($len >= $min_length) {
		print ">$id | $len\n" . $ctg_seq{$id} . "\n";
	}
}

exit(0);
