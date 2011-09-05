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
#  Module: Parser::FASTA::get_region_seq.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 fasta_file region_file\n";

my $fasta_file = shift or die $usage;
my $region_file = shift or die $usage;


open FASTAFILE, $fasta_file or die "Cannot open $fasta_file\n";
open REGIONFILE, $region_file or die "Cannot open $region_file\n";

my $seq = "";

### Get genome sequence
while( <FASTAFILE> ) {
	chomp;

	if (substr($_, 0, 1) eq ">") {
		# Do nothing, if just one chromosome
	}
	else {
		$seq .= $_;
	}
}

### Get regions
my $junk = <REGIONFILE>;

while( <REGIONFILE> ) {
	chomp;

	my @a = split(/\s+/, $_);
	
	my $reglen = $a[4] - $a[3] + 1;
	my $regseq = substr($seq, $a[3] - 1, $reglen);
	
	if($a[5] eq "-") {
		$regseq = revcomp($regseq);
	}

	print ">" . $a[0] . " | " . $a[1] . " | " . $a[2] . " | " . $a[3] . " | " . $a[4] . " | " . $a[5] . "\t";
	if   ($a[6] == 1) { print "tRNA\n"; }
	elsif($a[7] == 1) { print "rRNA\n"; }
	elsif($a[8] == 1) { print "sRNA\n"; }

	print "$regseq\n";
}

exit(0);

sub revcomp {
	my ($seq) = @_;

	my $newseq = reverse $seq;
	$newseq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

	return $newseq;
}

