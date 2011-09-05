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
#  Module: Analysis::CNV::parse_origin.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 originfile maplist\n";

my $file = shift or die $usage;
my $map = shift or die $usage;

###################################################################
# Read origins

my %O_POS = ();

open FILE, $file or die "Cannot open file\n";
while (<FILE>) {
	my @a = split " ";
	for (my $i = $a[1]; $i<=$a[2]; $i++) {
		my $id = $a[0] * 100000000 + $i;
		$O_POS{$id} = $a[0]."\t".$a[1]."\t".$a[2];
	}
}
close FILE;
print STDERR "Read origins\n";

####################################################################
# Read reads in origins

my %O_READ = ();

open MAP, $map or die $usage;
while (<MAP>) {
	my @a = split " ";
	if ($a[6] == 1 and ($a[9] == 4 or $a[9] == 7 or $a[9] == 10 or $a[9] == 13 or $a[9] == 16 or $a[9] == 19)) {
		my $id =  $a[0]*100000000 + $a[1];
		if (defined($O_POS{$id})) {
			$O_READ{$a[3]} = $O_POS{$id};
		}
	}
}
close MAP;
print STDERR "First pass\n";

####################################################################
# Read reads that are connected to origin reads

my %O_2_ANCHOR = ();

open MAP, $map or die $usage;
while (<MAP>) {
        my @a = split " ";
	if ($a[6] == 1 and ($a[9] == 4 or $a[9] == 7 or $a[9] == 10 or $a[9] == 13 or $a[9] == 16 or $a[9] == 19)) {
                my $id =  $a[0]*100000000 + $a[1];
                if (not defined($O_POS{$id}) and defined($O_READ{$a[3]})) {
                        $O_2_ANCHOR{$O_READ{$a[3]}} .= $a[0]."\t".$a[1]."\t".$a[4]."\n";
                }
        }
}	
close MAP;

print STDERR "Second pass\n";

####################################################################
# print

foreach my $origin (keys %O_2_ANCHOR) {
	print "#\t", $origin, "\n";
	print $O_2_ANCHOR{$origin};
}


