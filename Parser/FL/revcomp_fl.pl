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
#  Module: Parser::FL::revcomp_fl.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "reverse_fl.pl fl-file";
my $file = shift or die $usage;


open FILE, $file or die $usage;

while (<FILE>) {
	my @a = split " ", $_;
	print $a[0], "\t";
	print rev_comp($a[1]), "\t";
	for (my $i = 2; $i < @a; $i++) {
		my $s = reverse $a[$i];
		print $s;
		print "\t" if $i != @a -1;
	}
	print "\n";
}



sub rev_comp {
	my $seq = shift;
	my $rev = "";
	
	for (my $i = 0; $i<length($seq); $i++) {
		$rev .= rev_comp_base(substr($seq, $i, 1));
	}
	
	my $s = reverse $rev;

	return $s;
}

sub rev_comp_base {
        my ($base) = @_;

        $base =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

        return $base;
}




