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
#  Module: Analysis::SNP_by_Ref_qual::snp_by_ref_qual.pl
#  Purpose:
#  In:
#  Out:
#




open SNP, $ARGV[0];
open QUAL, $ARGV[1];

my $chr = 0;
my $pos = 1;
my %POS2QUAL = ();

QUALITY: while (my $line = <QUAL>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		# debug
		if ($line ne ">scaffold_1") {
			last QUALITY;
		}
		print STDERR $line, "\n";
		$pos = 1;	
		$chr++;
	}
	else {
		my @a = split " ", $line;
		for (my $i = 0; $i< @a; $i++) {
			$POS2QUAL{$chr."#".$pos} = $a[$i];
			$pos++;
		}
	}
}


close QUAL;


print STDERR "Let's go\n";

while (<SNP>) {
	my @a = split " ", $_;

	if (defined( $POS2QUAL{$a[0]."#".$a[1]})) {	
		print $a[0], "\t", $a[1], "\t", $POS2QUAL{$a[0]."#".$a[1]}, "\n";
	}
	else {
		die("Should not happen\n");
	}
}


