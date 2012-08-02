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
#  Module: Parser::FASTA::fasta_oneliner.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 min_length fastafile\n";

my $len = shift or die $usage;
my $in  = shift or die $usage;

open IN, $in or die $usage;

my $seq = "";
my $id  = "";
my %genome = ();

while (my $line = <IN>) {
	chomp($line);

	if (substr($line, 0, 1) eq ">") {

		if ($seq ne "") {
			$genome{$id} = $seq;
		}
		$seq = "";
		$id = substr($line, 1);
	}
	else {
		$seq .= $line;
	}
}

if ($seq ne "") {
	$genome{$id} = $seq;
}

### Find homo
foreach my $id_key ( sort keys %genome) {

	$seq=uc($genome{$id_key});

	while ($seq =~ /(A{$len,}|C{$len,}|G{$len,}|T{$len,})/g) {
		my $homopolymer_length = length($1);
		my $beg = length($`);
		my $end = $beg + $homopolymer_length - 1;
		print "$id_key\t$beg\t$end\t$1\n";
	}
}

