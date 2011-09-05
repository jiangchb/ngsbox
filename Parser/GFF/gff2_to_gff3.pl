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
#  Module: Parser::GFF::gff2_to_gff3.pl
#  Purpose:
#  In:
#  Out:
#


my $file = shift;

open IN, $file or die "Cannot open input file\n";
my $gene_id = "";

while(<IN>) {
	chomp;

	if( $_ =~ /gene/ ) {
		print "$_\n";
    		my @features = split(/\t/, $_);
		my ($gene_id, $junk) = split(/;/, $features[8]);
		$gene_id =~ s/ID=//;

		print "$features[0]\t$features[1]\tmRNA\t$features[3]\t$features[4]\t$features[5]\t" .
			"$features[6]\t$features[7]\tID=$gene_id.1;Parent=$gene_id\n";
	}
	elsif( $_ =~ /exon/ || $_ =~ /CDS/ ) {
		print "$_.1\n";
	}
	else { print STDERR "$_\n"; }

}

exit(0);
