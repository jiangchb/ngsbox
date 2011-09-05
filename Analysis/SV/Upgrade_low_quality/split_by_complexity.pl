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
#  Module: Analysis::SV::Upgrade_low_quality::split_by_complexity.pl
#  Purpose:
#  In:
#  Out:
#


use lib "$ENV{PGSP}/Analysis/SV";
use compare;

my $usage = "\nThis script splits the set of sv deletions into four disjunct set of SV deletions:
1. The overlapping ones (assumed to come from highly complex rearrangements)
2. The ones being completely included in others
3. Those common with the background
4. The very short ones
5. All others (assumned to be the \"simple\" deletions)

$0 LibID LibType MinSize SVFolder BackgroundSVFolder\n\n";

my $lib = shift or die $usage;
my $type = shift or die $usage;
my $size = shift or die $usage;
my $sv_folder = shift or die $usage;
my $bg_sv_folder = shift or die $usage;

my %SV = ();

###########################################################

print STDERR "Read background quality deletions...\n";

my $bgdelfile = $bg_sv_folder."/SV_deletion_high_quality_upgrade.lib".$lib.".".$type.".txt";
open FILE, $bgdelfile or die "Cannot open file: ".$bgdelfile."\n";

my @bgdel = ();

while (my $line = <FILE>) {
	my @a = split " ", $line;
	push @bgdel, $a[3]."#".$a[4]."#".$a[5];
}

close FILE;

print STDERR "...done\n";

###########################################################

print STDERR "Read quality deletions...\n";

my $delfile = $sv_folder."/SV_deletion_high_quality_upgrade.lib".$lib.".".$type.".txt";
open FILE, $delfile or die "Cannot open file: ".$delfile."\n";

my @dellines = ();

while (my $line = <FILE>) {
        push @dellines, $line;
}

close FILE;

print STDERR "...done\n";

###########################################################


print STDERR "Analyzing ", $sv_folder, "\n";

# Get min and max of insert dist
my ($min, $max) = get_min_max($sv_folder);
my $comp_obj = new compare();	

# Set files
my $complex_file = $sv_folder."/SV_deletion_high_quality_upgrade_complex.lib".$lib.".".$type.".txt";
my $background_file = $sv_folder."/SV_deletion_high_quality_upgrade_background.lib".$lib.".".$type.".txt";
my $short_file = $sv_folder."/SV_deletion_high_quality_upgrade_short.lib".$lib.".".$type.".txt";
my $simple_file = $sv_folder."/SV_deletion_high_quality_upgrade_simple.lib".$lib.".".$type.".txt";
my $included_file = $sv_folder."/SV_deletion_high_quality_upgrade_included.lib".$lib.".".$type.".txt";

open COMPLEX, ">".$complex_file or die "Cannot open file: ".$complex_file."\n";
open BG, ">".$background_file or die "Cannot open file: ".$background_file."\n";
open SHORT, ">".$short_file or die "Cannot open file: ".$short_file."\n";
open SIMPLE, ">".$simple_file or die "Cannot open file: ".$simple_file."\n";
open INCLUDED, ">".$included_file or die "Cannot open file: ".$simple_file."\n";


################################
## parse out overlapping
my @lines = ();
my $reg_chr = -1;
my $reg_end = -1;
my $ptr = 0;
foreach my $line (@dellines) {
	my @a = split " ", $line;

	# Set bg deletions to current
        while (is_smaller($bgdel[$ptr], $a[3], $a[4], $a[5]) == 1) {
		$ptr++;
		if (@bgdel < $ptr) {
                	last;
		}
	}

	# check for svs present in the background
	my $background = 0;
	my $num_overlapping = 0;
	while (overlaps($bgdel[$ptr+$num_overlapping], $a[3], $a[4], $a[5]) == 1) {
        	my @b = split "#", $bgdel[$ptr+$num_overlapping];
                if (compare->comp($min, $max, $a[3], $a[4], $a[5], $b[0], $b[1], $b[2]) == 1) {
                	$background = 1;
		}
                $num_overlapping++;
        }

	# print region 
	if ($background == 0) {
		if ((@lines > 0 and $a[4] > $reg_end) or ($a[3] != $reg_chr)) {
			print_region();
			@lines = ();
		}
	
		push @lines, $line;
		if ($a[5] > $reg_end || $a[3] != $reg_chr) {
			$reg_end = $a[5];
		}
		$reg_chr = $a[3];
	}
	else {
		print BG $line;
	}
}
print_region();


sub print_region {
	if (@lines == 0) {
		return;
	}
	elsif (@lines == 1) {
		my @a = split " ", $lines[0];
		if ($a[6] >= $size and $a[7] >= 10) {
			print SIMPLE $lines[0];
		}
		else {
			print SHORT $lines[0];
		}
	}
	else {
		# set beg and end
		my $beg = 0; my $end = 0;
		for (my $i = 0; $i < @lines; $i++) {
			my @a = split " ", $lines[$i];
			if ($a[4] < $beg || $beg == 0) {
				$beg = $a[4];
			} 
			if ($a[5] > $end || $end == 0) {
				$end = $a[5];	
			}
		}
		# one overlaps all ?
		for (my $i = 0; $i < @lines; $i++) {
			my @a = split " ", $lines[$i];
                        if ($a[4] == $beg and $a[5] == $end) {
				for (my $j=0; $j < $i; $j++) {
					print INCLUDED $lines[$j];
				}
				if ($a[6] >= $size and $a[7] >= 10) {
					print SIMPLE $lines[$i];
				}
				else {
					print SHORT $lines[$i];
				}
				for (my $j=$i+1; $j < @lines; $j++) {
					print INCLUDED $lines[$j];
				}
				return;
			}
		}
		# else: it is a complex region
		for (my $i = 0; $i < @lines; $i++) {
			print COMPLEX $lines[$i];
		}
	}
}

sub overlaps {
        my ($pat, $chr, $pos, $end) = @_;
        my @a = split "#", $pat;

        if ($a[0] != $chr) {
                return 0;
        }
        elsif (($a[1] >= $pos and $a[1] <= $end) or ($pos >= $a[1] and $pos <= $a[2])) {
                return 1;
        }

        return 0;
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

sub is_smaller {
        my ($pat, $chr, $pos, $end) = @_;
        my @a = split "#", $pat;

        if ($a[0] < $chr) {
                return 1;
        }
        elsif ($a[0] == $chr and $a[2] < $pos) {
                return 1;
        }

        return 0;
}


