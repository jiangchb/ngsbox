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
#  Module: Parser::ML::separate_mappings_by_rep_state.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "perl separate_mappings_by_rep_state.pl map.list\n";
my $file = shift or die $usage;
open FILE, $file;

my %FIRST = ();
my %SECOND = ();

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $rep_state = "";
	if ($a[6] == 1) {
		$rep_state = 1;
	}
	else {
		$rep_state = 0;
	}

	if ($a[9] == 1 or $a[9] == 4 or $a[9] == 5 or $a[9] == 6 or $a[9] == 9 or $a[9] == 10 or$a[9] == 11) {
		$FIRST{$a[3]} = $rep_state;
	}
	else {
		$SECOND{$a[3]} = $rep_state;
	}
}

close FILE;

open FILE, $file;
open REPREP, ">".$file.".rep-rep";
open REPUNIQ, ">".$file.".rep-uniq";
open UNIQUNIQ, ">".$file.".uniq-uniq";

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (defined($FIRST{$a[3]}) and defined($SECOND{$a[3]})) {
		if ($FIRST{$a[3]} == 0 and $SECOND{$a[3]} == 0) {
			print REPREP $line;
		}
		if ($FIRST{$a[3]} == 1 and $SECOND{$a[3]} == 0) {
                        print REPUNIQ $line;
                }
		if ($FIRST{$a[3]} == 0 and $SECOND{$a[3]} == 1) {
                        print REPUNIQ $line;
                }
		if ($FIRST{$a[3]} == 1 and $SECOND{$a[3]} == 1) {
                        print UNIQUNIQ $line;
                }
	}
}




