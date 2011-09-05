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
#  Module: Analysis::TE::check_genes_4_presence.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 SVdelfile reffile varfile tairfile sample\n";

my $svdel = shift or die $usage;
my $ref = shift or die $usage;
my $var = shift or die $usage;
my $tair = shift or die $usage;
my $sample = shift or die $usage;

######################################################################
# Missing accession data
######################################################################

my %DEL=();

open DEL, $svdel;
while (<DEL>) {
	my @a = split " ";
	for (my $i = $a[4]; $i<= $a[5]; $i++) {
		$DEL{$a[3]."#".$i} = 1;
	}
}
close DEL;

print STDERR "LOADED ALL DELETIONS\n";


######################################################################
#  Variation poses something like presence
######################################################################

my %PRESENCE = ();
open REF, $ref;
while (<REF>) {
        my @a = split " ";
       	$PRESENCE{$a[1]."#".$a[2]} = 1;
}
close REF;

open VAR, $var;
while (<VAR>) {
        my @a = split " ";
        $PRESENCE{$a[1]."#".$a[2]} = 1;
}
close VAR;


######################################################################

open TAIR, $tair;

while (<TAIR>) {
	chomp();
	my @a = split " ";

	my $count_sv = 0;
	my $count_presence = 0;

	for (my $i = $a[3]; $i<= $a[4]; $i++) {
		if (defined($DEL{$a[2]."#".$i})) {
			$count_sv++;
		}
		if (defined($PRESENCE{$a[2]."#".$i})) {
                        $count_presence++;
                }
	}

	print $sample, "\t", $a[0], "\t", $a[1], "\t", $a[4]-$a[3]+1, "\t", $count_presence, "\t", $count_sv, "\n";
}


