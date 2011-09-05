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
#  Module: Parser::ML::select_read_from_region.pl
#  Purpose:
#  In:
#  Out:
#


my $map_list = shift;
my $chr      = shift;
my $start    = shift;
my $end      = shift;
my $type     = shift; # fasta oder fastq
 

### Get reads from target region
open FILE, $map_list or die "cannnot open $map_list\n";
my $written = 0;
my %READS = ();
while (<FILE>) {
	my @a = split " ", $_;

	if( $a[0] == $chr && $a[1] >= $start && $a[1] <= $end) {
		if (not defined($READS{$a[3].".".$a[9]})) {
			if ($type eq "fasta") {
				print ">", $a[3], "\n", get_seq($a[2]), "\n";
			}
			elsif ($type eq "fastq") {
				print "@", $a[3], "\n", get_seq($a[2]), "\n+\n", $a[11], "\n";
			}
			$READS{$a[3].".".$a[9]} = 1;
		}
		$written = 1;
	}
	else {
		if ($written == 1) {
			exit(0);
		}
	}
}


exit(0);

sub get_seq {
	my ($align) = @_;
	my $seq = "";
        for (my $i = 0; $i<length($align); $i++) {
		if (substr($align, $i, 1) eq "[") {
			if (substr($align, $i+2, 1) ne "-") {
	                	$seq .= substr($align, $i+2, 1);
			}
                	$i+=3;
                }
                else {
                	$seq .= substr($align, $i, 1);
                }
        }

        return $seq;
}

