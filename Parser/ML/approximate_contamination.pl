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
#  Module: Parser::ML::approximate_contamination.pl
#  Purpose:
#  In:
#  Out:
#


my %R = ();

while (<STDIN>) {
	my ($chr, $pos, $alg, $id) = split " ";	
	$R{$id} .= $chr."#";
}

my $c = 0;
my $c_g = 0;
my $c_o = 0;

foreach my $val (values %R) {
	my @chrs = split "#", $val;
	my $o_flag = 0;
	my $g_flag = 0;
	foreach my $i (@chrs) {
		if ($i==702 || $i==696 || $i==697 || $i==698 || $i==699 || $i==700 || $i==701) {
			$o_flag = 1;
		} else {
			$g_flag = 1;
		}
	}	
	# set counts
	if ($o_flag == 0) {
		$c_g++;
	} elsif ($g_flag == 0) {
		$c_o++;
	}
	$c++;
}

print STDOUT "all:\t\t", $c, "\n";
print STDOUT "ambi:\t\t", $c-($c_g+$c_o), "\n";
print STDOUT "genomic:\t", $c_g, "\n";
print STDOUT "organelles:\t", $c_o, "\n";



