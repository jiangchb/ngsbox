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
#  Module: Analysis::Eighties::CNV::anno2shore80.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 NBSLRR TAIR8_genes_trasposons\n";

my $base = "/ebio/abt6_projects/backup/data/solexa_analysis/ATH/Genome/Col-0-Geneva/AlignmentFolder/CNV_per_segment/";

my $nbslrr = shift or die $usage;
my $tair8 = shift or die $usage;

#########################################################################################
my %NBS = ();
open FILE, $nbslrr or die "Cannot open file\n";
while (my $line = <FILE>) {
	chomp($line);
	$NBS{$line} = 1;
}
close FILE;

#########################################################################################
my %POS2ANNO = ();
open FILE, $tair8 or die "Cannot open file\n";
while (my $line = <FILE>) {
	if (substr($line, 0, 1) ne "#") {
		my @a = split " ", $line;
		$POS2ANNO{$a[0]."#".$a[1]} = $a[4];
	}
}
close FILE;

#########################################################################################

my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE173","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");

my %ID = ();

my @ECO_VALUE = {};
my %COL_VALUE = {};
my @ECO_COUNT = ();
my $COL_COUNT = 0;

for (my $eco = 0; $eco < @ECOTYPES; $eco++) {
	$COL_COUNT = 0;
	my $eco_count = 0;
	my %eco_value = {};

	my $file = $base."/shore_count.".$ECOTYPES[$eco].".txt";

	open FILE, $file or die "Cannot open file: ".$file."\n";

	while (my $line = <FILE>) {
		if (not(substr($line, 0, 1) eq "#") and length($line) > 5)  {
			my @a = split " ", $line;

	                chomp($line);
		
			$COL_COUNT += $a[4];
			$eco_count += $a[5];
			
			my $id = ($a[0]*100000000)+$a[1];

			$COL_VALUE{$id} = $a[4];
			$eco_value{$id} = $a[5];
			$ID{$id} = $a[1]+$a[2];
		}	
	}
	
	$ECO_VALUE[$eco] = \%eco_value;
	$ECO_COUNT[$eco] = $eco_count;

	print STDERR $ECOTYPES[$eco], "\t", $eco_count, "\t", (keys %eco_value)+0, "\n";

	close FILE;
}

print @ECO_VALUE+0, "\n";
print @ECO_COUNT+0, "\n";

#for (my $i = 0; $i < @ECO_VALUE; $i++) { print $i, "\t"; print ((keys %{$ECO_VALUE[$i]})+0); print "\n"; }

#########################################################################################
## Output 

my @OUT_LINES = ();

# init lines with lines pos and col-0 values
foreach my $id (sort {$a <=> $b}  keys %ID) {
	my $chr = int($id/100000000);
	my $pos = $id%100000000;
	my $end = $ID{$id};
	my $col_val = $COL_VALUE{$id};
	my $col_norm = $col_val / ($COL_COUNT / 1000000);

	my $anno = "\\N";
	my $func = "\\N";

	if (defined($POS2ANNO{$chr."#".$pos})) {
		$anno = $POS2ANNO{$chr."#".$pos};
		if (defined($NBS{$POS2ANNO{$chr."#".$pos}})) {
			$func = "NBS_LRR_active"; 
		}
	}

	my $line = $anno."\t".$func."\t".$chr."\t".$pos."\t".$end."\t".$col_val."\t".$col_norm;

	push @OUT_LINES, $line;
}

print STDERR "Initialized lines\n";

# for each ecotype parse over lines and add repective value
for (my $eco = 0; $eco < @ECOTYPES; $eco++) {
	my %eco_value = %{$ECO_VALUE[$eco]};
	my $eco_count = $ECO_COUNT[$eco];
	print STDERR $ECOTYPES[$eco], "\t", $eco_count, "\t", (keys %eco_value)+0, "\n";
	for (my $line = 0; $line < @OUT_LINES; $line++) {
		my @a = split " ", $OUT_LINES[$line];
		my $id = (100000000*$a[2])+$a[3];
		my $norm = $eco_value{$id} / ($eco_count / 1000000);
		my $frac = 0;
		my $col_val = $COL_VALUE{$id};
		my $col_norm = $col_val / ($COL_COUNT / 1000000);

		if ($norm < $col_norm) {
			if ($norm != 0) {
				$frac = ($col_norm / $norm) * (-1);
			}
			else {
				$frac = $col_norm * (-1);
			}
			$frac = -7 if $frac < -7;
			$frac = -1 if $frac > -1;
		}
		else {
			if ($col_norm != 0) {
				$frac = ($norm / $col_norm);
			}
			else {
				$frac = $norm;
			}
                        $frac = 7 if $frac > 7;
			$frac = 1 if $frac < 1;
		}

		$OUT_LINES[$line] .= "\t".$eco_value{$id}."\t".$norm."\t".$frac;
	}
}	


# print out
open OUTFILE, ">".$base."/shore_count.ALL_80.Bak-7.notcomplete.txt";
for (my $line = 0; $line < @OUT_LINES; $line++) {
	print OUTFILE $OUT_LINES[$line], "\n";
}	

close FILE;



