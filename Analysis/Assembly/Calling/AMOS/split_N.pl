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
#  Module: Analysis::Assembly::Calling::AMOS::split_N.pl
#  Purpose:
#  In:
#  Out:
#



### User params
my $contig_fasta = shift;
open CONTIG, $contig_fasta or die "Cannot open $contig_fasta\n";

my $counter = 1;
while(<CONTIG>) {
	chomp($_);

	if (substr($_, 0, 1) eq ">") {
		my $id = substr($_, 1);
		my $seq = <CONTIG>;
		chomp($seq);

		if( ($seq =~ /NNN/) && ($id =~ /velvet/) ) {
			my @sub_seqs = split(/N+/, $seq);
			foreach my $sub_seq (@sub_seqs) {
				print ">$counter$id\n$sub_seq\n";
				$counter++;
			}
		}
		else {
			print ">$id\n$seq\n";
		}
	}
}
close CONTIG;
