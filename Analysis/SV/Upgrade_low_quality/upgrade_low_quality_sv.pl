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
#  Module: Analysis::SV::Upgrade_low_quality::upgrade_low_quality_sv.pl
#  Purpose:
#  In:
#  Out:
#


use lib "$ENV{PGSP}/Analysis/SV";
use compare;

my $usage = "This script upgrades low quality sv deletions if high quality deletions
exist in a different strain, and if the deletion to be upgraded has at leat
MinUncovered percent (as double value) positions.
Output a combined set of high quality and upgraded deletions.
$0 LibID LibType MinUncovered SVFolder SVFolder [SVFolder [SVFolder [...]]]\n\n";

my $lib = shift or die $usage;
my $type = shift or die $usage;
my $uncovered = shift or die $usage;

my %SV = ();

print STDERR "Read high quality deletions...\n";

foreach my $sv_folder (@ARGV) {
	my $file = $sv_folder."/SV_deletion_high_quality.lib".$lib.".".$type.".txt";
	open FILE, $file or die "Cannot open file: ".$file."\n";

	my $acc = "";
	my @del = ();

	while (my $line = <FILE>) {
		my @a = split " ", $line;
		$acc = $a[0];
		push @del, $a[3]."#".$a[4]."#".$a[5];
	}
	$SV{$acc} = \@del;

	close FILE;
}

print STDERR "...done\n";

###########################################################

my %SV_ptr = ();


foreach my $sv_folder (@ARGV) {

	print STDERR "Analyzing ", $sv_folder, "\n";

	# Get min and max of insert dist
	my ($min, $max) = get_min_max($sv_folder);
	my $comp_obj = new compare();	

	# Set files
	my $file = $sv_folder."/SV_deletion_low_quality.lib".$lib.".".$type.".txt";
	my $high = $sv_folder."/SV_deletion_high_quality.lib".$lib.".".$type.".txt";
        open FILE, $file or die "Cannot open file: ".$file."\n";
	my $out = $sv_folder."/SV_deletion_high_quality_upgrade.lib".$lib.".".$type.".tmp";
	open UPGRADE, ">".$out or die "Cannot open file: ".$out."\n";

	set_ptr();

	while (my $line = <FILE>) {
		my @a = split " ", $line;

                my $acc = $a[0];
		my $chr = $a[3];
		my $pos = $a[4];
		my $end = $a[5];

		SV: foreach my $cmp_acc (keys %SV) {
			if ($acc ne $cmp_acc) {
				# skip all which are upstream and thus cannot overlap anymore
#print STDERR ${$SV{$cmp_acc}}[$SV_ptr{$cmp_acc}], "\t", $chr, "\t", $pos, "\t", $end, "\n";
				while (is_smaller(${$SV{$cmp_acc}}[$SV_ptr{$cmp_acc}], $chr, $pos, $end) == 1) {
					$SV_ptr{$cmp_acc}++;
					if (@{$SV{$cmp_acc}} < $SV_ptr{$cmp_acc}) {
						last;
					}
				}
				# Compare with all overlapping SVs
				my $num_overlapping = 0;
				while (overlaps(${$SV{$cmp_acc}}[$SV_ptr{$cmp_acc}], $chr, $pos, $end) == 1) {
					my @b = split "#", ${$SV{$cmp_acc}}[$SV_ptr{$cmp_acc}];		
					if (compare->comp($min, $max, $chr, $pos, $end, $b[0], $b[1], $b[2]) == 1) {
						if ($a[8]/$a[7] >= $uncovered) {
							print UPGRADE $line;
							last SV;
						}
					}
					$SV_ptr{$cmp_acc}++;
					$num_overlapping++;
				}
				$SV_ptr{$cmp_acc}-=$num_overlapping; # set back that next SV has also the chance to be compared with these guys
			}
		}
	}

	close FILE;
	close UPGRADE;

	####### 
	my $cmd = "sort -m -k4n -k5n -k6n $high $out > ".$sv_folder."/SV_deletion_high_quality_upgrade.lib".$lib.".".$type.".txt";
	system($cmd);
	system("rm $out")
	

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

sub set_ptr {
	foreach my $acc (keys %SV) {
		$SV_ptr{$acc} = 0;
	}
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


