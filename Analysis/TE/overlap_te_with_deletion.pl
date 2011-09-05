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
#  Module: Analysis::TE::overlap_te_with_deletion.pl
#  Purpose:
#  In:
#  Out:
#


my $delfile = "/ebio/abt6_projects/backup/data/project_analysis/TE_effects_on_expression/DataAquisition/Deleted_TEs/Bur-0.deletion.tab" or die "Cannot open file\n";
my $delfile2 = "/ebio/abt6_projects/backup/data/project_analysis/TE_effects_on_expression/DataAquisition/Deleted_TEs/C24.deletion.tab" or die "Cannot open file\n";
my $tefile = "/ebio/abt6_projects/backup/data/project_analysis/TE_effects_on_expression/DataAquisition/Deleted_TEs/te.tab" or die "Cannot open file\n";

my %DEL=();
open DEL, $delfile or die "Cannot open file\n";
while (<DEL>) {
	my @a = split " ";
	for (my $i = $a[4]; $i<= $a[5]; $i++) {
		$DEL{$a[3]."#".$i} = 1;
	}
}
close DEL;


my %DEL2=();
open DEL2, $delfile2 or die "Cannot open file\n";
while (<DEL2>) {
        my @a = split " ";
        for (my $i = $a[4]; $i<= $a[5]; $i++) {
                $DEL2{$a[3]."#".$i} = 1;
        }
}
close DEL2;


print STDERR "LOADED ALL DELETIONS\n";

open TE, $tefile or die "Cannot open file\n";

while (<TE>) {
	chomp();
	my @a = split " ";
	my $count = 0;
	my $count2 = 0;
	my $length = 0;
	for (my $i = $a[2]; $i<= $a[3]; $i++) {
		$length++;
		if (defined($DEL{$a[1]."#".$i})) {
			$count++;
		}
		if (defined($DEL2{$a[1]."#".$i})) {
                        $count2++;
                }
	}
	print $_, "\t", $count, "\t", $count2, "\n";
}


