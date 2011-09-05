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
#  Module: Parser::FASTQ::rev_comp.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 fastqfile\n";
my $file = shift or die $usage;

open FILE, $file or die $usage;

while (my $line = <FILE>) {
	print $line;
	if (substr($line, 0, 1) eq "@" or substr($line, 0, 1) eq "+") { 
		my $lline = <FILE>;
		chomp($lline);
		if (substr($line, 0, 1) eq "@") {
			print rev_comp($lline), "\n";
		}
		else {
			$lline = reverse($lline);
			print $lline, "\n";
		}
	}
}

sub rev_comp {
        my ($seq) = @_;

        $seq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;
	$seq = reverse($seq);

        return $seq;
}




