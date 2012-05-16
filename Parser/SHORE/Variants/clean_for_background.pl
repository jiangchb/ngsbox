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
#  Module: Parser::SHORE::Variants::clean_for_background.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 background,background,... variants\n\n";

my $bg = shift or die $usage;
my $var = shift or die $usage;

my %BG = ();
my @VAR = ();

my @files = split ",", $bg;

foreach my $file (@files) {

	open FILE, $file or die $usage;

	while (my $line = <FILE>) {
		my @a = split " ", $line;
		$BG{$a[1]."#".$a[2]} = 1;
	}

	close FILE;
}


open FILE, $var or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (not defined($BG{$a[1]."#".$a[2]})) {
	#	my $id = ($a[1] * 100000000) + $a[2];
		push @VAR, $line;
	}

}

foreach my $key (@VAR) {
	print $key;
}



