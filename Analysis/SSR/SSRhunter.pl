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
#  Module: Analysis::SSR::SSRhunter.pl
#  Purpose:
#  In:
#  Out:
#



use lib "$ENV{PGSP}/Prediction/MA/";

use SSR;

my @files = ("mono.txt", "di.txt", "tri.txt", "tetra.txt", "penta.txt", "hexa.txt");
my $max_length = 28;
my $vicinity_length = 3;

my $mapfile = shift or die "Set map.list file\n";
my $maline = shift or die "Set ma line\n";
my $chromosome = shift or die "Set chromosome\n";

my %SSRS = ();

foreach my $file (@files) {
	open FILE, $file or die "Cannot open file\n";

	while (<FILE>) {
		my ($chr, $pos, $end, $len, $ins, $seq) = split " ";

		if ($len <= $max_length && $chr == $chromosome) {
			my $id = (100000000 * $chr) + $pos;
			$SSRS{$id} = new SSR();

			$SSRS{$id}->set($chr, $pos, $end, $len, $ins, $seq);	

		}
	}

	close FILE;
}


open MAP, $mapfile or die "Cannot open mapfile\n";

my @ssr_ids = sort {$a <=> $b} keys %SSRS;

print STDERR "Found: ", @ssr_ids+0, "\n";;

my $ssr_id = 0;
my $ssr_chr = $SSRS{$ssr_ids[$ssr_id]}{chromosome};
my $ssr_pos = $SSRS{$ssr_ids[$ssr_id]}{position};

my $count = 0;
MAPPINGS: while (<MAP>) {
	my ($chr, $pos, $alignment, $readid, $strand, $mm, $hits, $readlength, $pe, $qual1, $qual2, $qual3) = split " ";
	
	$count++;
	
	while ($chr > $ssr_chr || $pos >= $ssr_pos) {
		$ssr_id++;
		if ($ssr_id >= @ssr_ids) {
			last MAPPINGS;
		}
		$ssr_chr = $SSRS{$ssr_ids[$ssr_id]}{chromosome};
		$ssr_pos = $SSRS{$ssr_ids[$ssr_id]}{position};
	}

	my $flanking_length = $readlength-$SSRS{$ssr_ids[$ssr_id]}{len};

	if ($pos >= $ssr_pos-($flanking_length)+3  && $pos < $ssr_pos-3) {
		if ($alignment !~ /-/) {
			$SSRS{$ssr_ids[$ssr_id]}{spanning_ungapped}++;
		}
		else {
			$SSRS{$ssr_ids[$ssr_id]}{spanning_gapped}++;
			$SSRS{$ssr_ids[$ssr_id]}{variation} .= 	get_variation($alignment).","
		}
	}
}

for (my $i = 0; $i < @ssr_ids; $i++) {
	#if ($SSRS{$ssr_ids[$i]}{spanning} < 2) {
		$SSRS{$ssr_ids[$i]}{variation} = "-" if $SSRS{$ssr_ids[$i]}{variation} eq "";
		print $maline, "\t", $SSRS{$ssr_ids[$i]}{chromosome}, "\t", $SSRS{$ssr_ids[$i]}{position}, "\t", $SSRS{$ssr_ids[$i]}{end}, "\t", $SSRS{$ssr_ids[$i]}{len}, "\t", $SSRS{$ssr_ids[$i]}{instances}, "\t", $SSRS{$ssr_ids[$i]}{sequence}, "\t", $SSRS{$ssr_ids[$i]}{spanning_ungapped}, "\t", $SSRS{$ssr_ids[$i]}{spanning_gapped}, "\t", $SSRS{$ssr_ids[$i]}{variation}, "\n";
	#}
}


sub get_variation {
	my $alignment = $_;

	my $variation = "";
	my $flag = 0;
	for (my $i = 0; $i<length($alignment); $i++) {
		$flag = 0;
		while (substr($alignment, $i, 1) eq "[") {
			if (substr($alignment, $i+1, 1) eq "-" || substr($alignment, $i+2, 1) eq "-") {
				$flag = 1;
				$i++;
				if (substr($alignment, $i, 1) eq "-") {
					$i++;
					$variation .= substr($alignment, $i, 1);
				}
				else {
					$variation .= lc(substr($alignment, $i, 1));
					$i++;
				}
				$i+=2;
			}
			else {
				$i+=4;
			}
		}
		$variation .= "#" if $flag == 1;
	}

	
	return $variation;
}

