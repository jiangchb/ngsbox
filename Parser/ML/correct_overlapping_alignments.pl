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
#  Module: Parser::ML::correct_overlapping_alignments.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 map.list\n";
open FILE, shift or die $usage;
open OUT, ">map.list.corrected";

my %LINES = ();

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $pe_flag = $a[9];
	my $id = $a[3];
	my $begin_r2 = $a[1];
	if (is_happy($pe_flag)) {
		if (defined($LINES{$id})) {
			my $end1 = get_end($LINES{$id});
		
			if ($begin_r2 > $end1) {
				print OUT $LINES{$id};
				print OUT $line;
			}
			else {
				my $line2 = correct4overlap($end1, $LINES{$id}, $line); 
				if ($line2 ne "") {
					print OUT $LINES{$id};
					print OUT $line2;
				}
				else {
					print OUT $LINES{$id};
				}
			}

			delete $LINES{$id};
		}
		else {
			# store read
			$LINES{$id} = $line;
		}
	}	
	else {
		print OUT $line;
	}
}

foreach my $ll (values %LINES) {
	print OUT $ll;
}

close FILE;
close OUT;

system("sort -k1n -k2n map.list.corrected > map.list.corrected.sorted");
system("rm map.list.corrected");

sub correct4overlap {
	my ($end_read1, $line) = @_;

	my @a = split " ", $line;

	my $begin = $a[1];
	my $read_bases_skipped = 0;
	my $align = $a[2];
	my $i = 0;
	
	# set begin to the genomic position where the two alignments do not overlap anymore
	# and i to the number of characters that are skipped in the second alignment
	IT: while ($i <= length($align)) {
		if ($begin > $end_read1) {
			last IT;
		}

		if (substr($align, $i, 1) eq "[") {
			if (substr($align, $i+1, 1) ne "-") { # if ref eq "-" than the ref pos is not shifted forward
				$begin++;
			}
			if (substr($align, $i+2, 1) ne "-") {
				$read_bases_skipped++;
			}
			$i+=4;
		}
		else {
			$begin++;
			$read_bases_skipped++;
			$i+=1;
		}
	
	}

	# set new alignment and parse for new alignment specific values
	my $new_alignment = substr($align, $i, length($align)-$i);
	my $new_mm_count = 0;

	for (my $u = 0; $u < length($new_alignment); $u++) {
		if (substr($a[2], $u, 1) eq "[") {
                	$new_mm_count++;
                        $u+=3;
                }
        }
	

	my $new_line = "";
	$new_line = $a[0]."\t".$begin."\t";
	$new_line .= $new_alignment."\t";
	$new_line .= $a[3]."\t".$a[4]."\t".$new_mm_count."\t".$a[6]."\t";
	$new_line .= ($a[7]-$read_bases_skipped)."\t".$read_bases_skipped."\t".$a[9]."\t";

	if ($a[4] eq "D") {
		$new_line .= substr($a[10], $read_bases_skipped, length($a[10])-$read_bases_skipped)."\t".substr($a[11], $read_bases_skipped, length($a[11])-$read_bases_skipped)."\t".substr($a[12], $read_bases_skipped, length($a[12])-$read_bases_skipped)."\n";
	}
	else {
		$new_line .= substr($a[10], 0, length($a[10])-$read_bases_skipped)."\t".substr($a[11], 0, length($a[11])-$read_bases_skipped)."\t".substr($a[12], 0, length($a[12])-$read_bases_skipped)."\n";
	}


	# Alignment that are too short are skipped
	if (length($new_alignment) > 1 and ($a[7]-$read_bases_skipped) > 10) {
		return $new_line;
	}
	else {
		return "";
	}

}

sub get_end {
	my ($line) = @_;
	my @a = split " ", $line;
	my $refpos = 0;
	for (my $i = 0; $i < length($a[2]); $i++) {
		my $nuc = substr($a[2], $i, 1);
		if ($nuc eq "[") {
			if (substr($a[2], $i+1, 1) ne "-") {
				$refpos++;
			}
			$i+=3;
		}
		else {
			$refpos++;
		}
	}
	return $a[1]+$refpos-1;
}

sub is_happy {
	my ($flag) = @_;
	if ($flag == 3 or $flag == 6 or $flag == 9 or $flag == 12) {
		return 1;
	}
	return 0;
}


