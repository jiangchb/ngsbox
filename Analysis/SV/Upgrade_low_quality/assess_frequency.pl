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
#  Module: Analysis::SV::Upgrade_low_quality::assess_frequency.pl
#  Purpose:
#  In:
#  Out:
#


use lib "$ENV{PGSP}/Analysis/SV";
use compare;

my $usage = "This script assesses the frequency of deletions in a set
of SV deletion prediction. It uses either the simple or the complex
deletions. (First; set Complex to 0, or latter set Complex = 1
$0 LibID LibType Complex SVFolder SVFolder [SVFolder [SVFolder [...]]]\n\n";

my $lib = shift or die $usage;
my $type = shift or die $usage;
my $complex = shift;

if ($complex ne "1" and $complex ne "0") { die "Complex needs to be either 1 or 0."; }

my @SV_CHR = ();
my @SV_POS = ();
my @SV_END = ();
my @SV_ID = ();

my $ID = 0;

print STDERR "Read high quality deletions...\n";

foreach my $sv_folder (@ARGV) {
	my $file = $sv_folder."/SV_deletion_high_quality_upgrade_simple.lib".$lib.".".$type.".txt";
	if ($complex == 1) {  $file = $sv_folder."/SV_deletion_high_quality_upgrade_complex.lib".$lib.".".$type.".txt"; }

	open FILE, $file or die "Cannot open file: ", $file, "\n";
	
	my $output_folder = $sv_folder."/SV_deletion_high_quality_upgrade_simple_frequency_".(@ARGV+0);
	if ($complex == 1) {  $output_folder = $sv_folder."/SV_deletion_high_quality_upgrade_complex_frequency_".(@ARGV+0); }

	if (not -e $output_folder) {
		mkdir($output_folder);
	}

	my $outfile = ">".$output_folder."/SV_deletion_high_quality_upgrade_simple_id.lib".$lib.".".$type.".txt";
	if ($complex == 1) {  $outfile = ">".$output_folder."/SV_deletion_high_quality_upgrade_comlex_id.lib".$lib.".".$type.".txt"; }

	open OUT, $outfile or die "Cannot open file: ", $outfile, "\n";

	print STDERR "Analyzing ", $sv_folder, "\n";

	# Get min and max of insert dist
	my ($min, $max) = get_min_max($sv_folder);
	my $comp_obj = new compare();	

	while (my $line = <FILE>) {
		my @a = split " ", $line;

                my $acc = $a[0];
		my $chr = $a[3];
		my $pos = $a[4];
		my $end = $a[5];

		my $added = 0;

		for (my $i = 0; $i < @SV_CHR; $i++) {

			my $cmp_chr = $SV_CHR[$i];
			my $cmp_pos = $SV_POS[$i];
			my $cmp_end = $SV_END[$i];
			my $cmp_id = $SV_ID[$i];

			if ($cmp_chr == $chr and (($cmp_pos >= $pos and $cmp_pos <= $end) or ($pos >= $cmp_pos and $pos <= $cmp_end))) {
				# Compare with all overlapping SVs
				if (compare->comp($min, $max, $chr, $pos, $end, $cmp_chr, $cmp_pos, $cmp_end) == 1) {
				
					push @SV_CHR, $chr;
	                	        push @SV_POS, $pos;
        		                push @SV_END, $end;
					push @SV_ID, $cmp_id;

					$added = $cmp_id;
					last;
				}
			}
			
		}

		if ($added == 0) {
			$ID++;
			$added = $ID;
			push @SV_CHR, $chr;
			push @SV_POS, $pos;
			push @SV_END, $end;
			push @SV_ID, $ID;
		}

		chomp($line);
		print OUT $line, "\t", $added, "\n";

	}

	close FILE;
	close OUT;


}


sub get_min_max {
	my ($folder) = @_;
        my $tmpfile = "$folder/insert.tmp";
        system("~/shore/startshore.sh struct -I $folder/insert_dist.txt > $tmpfile");
        open TMP, $tmpfile;
        my $l = <TMP>;
        my @a = split " ", $l;
	return ($a[2], $a[3]);
}




