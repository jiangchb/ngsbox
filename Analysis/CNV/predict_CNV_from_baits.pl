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
#  Module: Analysis::CNV::predict_CNV_from_baits.pl
#  Purpose:
#  In:
#  Out:
#

use Statistics::Lite qw(:all);
use Getopt::Long;


### Command line parameters
# files
my $file_info;
my $file_classification;
my $file_uncertainty;
my $file_means;
my $file_proportions;

# parameters
my $foldchange = 1.3;
my $chr;

my %CMD;
GetCom();


#print "\n\n$file_info\n$file_classification\n$file_uncertainty\n$file_means\n$file_proportions\n$foldchange\n$chr\n\n";

### Parse info
my $i = 0;
my @info  = ();
open INFO, $file_info or die "Cannot open info file\n";
<INFO>;
while( <INFO> ) {
	chomp;
	my @a = split("\t", $_);
	$info[$i] = \@a;
	$i++;
}
close INFO or die;


### Parse classification (sample vs. bait)
$i = 0;
my @classification = ();
open CLASSI, $file_classification or die "Cannot open classification file\n";
<CLASSI>;
while( <CLASSI> ) {
	chomp;
	my @a = split("\t", $_);
	$classification[$i] = \@a;
	$i++;
}
close CLASSI or die;


### Parse uncertainty (sample vs. bait)
# $i = 0;
# my @uncertainty = ();
# open UNCERT, $file_uncertainty or die "Cannot open uncertainty file\n";
# <UNCERT>;
### -> use push to invert this matrix
# while( <UNCERT> ) {
# 	chomp;
# 	my @a = split("\t", $_);
# 	$uncertainty[$i] = \@a;
# 	$i++;
# }
# close UNCERT or die;


### Parse uncertainty (bait vs. sample)
my @uncertainty = ();
open UNCERT, $file_uncertainty or die "Cannot open uncertainty file\n";
<UNCERT>;
while( <UNCERT> ) {
	chomp;
	my @a = split("\t", $_);
	for(my $x = 0; $x <= $#a; $x++) {
		push( @{$uncertainty[$x]}, $a[$x] );
	}
}
close UNCERT or die;

#print join("\n", @{$uncertainty[8]});
#exit(0);

### Parse mean (bait vs. group-mean-cov)
$i = 0;
my @mean = ();
open MEAN, $file_means or die "Cannot open mean file\n";
<MEAN>;
while( <MEAN> ) {
	chomp;
	my @a = split("\t", $_);
	$mean[$i] = \@a;
	$i++;
}
close MEAN or die;


### Parse proportions (sample vs. bait)
$i = 0;
my @proportions = ();
open PROP, $file_proportions or die "Cannot open proportions file\n";
<PROP>;
while( <PROP> ) {
	chomp;
	my @a = split("\t", $_);
	$proportions[$i] = \@a;
	$i++;
}
close PROP or die;


### Parse matrix to predict consecutive CNV paterns
my %cluster = (); # store consecutive baits that fulfill the CNV criteria and the pattern similarity of >=95%
my %cluster_fc = ();
my %cluster_uc = ();
my $cluster_start_id = -10;
my $lastid = -10;

for(my $n = 0; $n <= $#mean; $n++) {

	### If at least two CNV allele groups exist
	if( $mean[$n][1] ne "NA" ) {

		### Sum up rpkm of all groups
		my $sum_rpkm = 0;
		for(my $g = 0; $g < 9; $g++) {
			if($mean[$n][$g] ne "NA") {
				$sum_rpkm += $mean[$n][$g];
			}
		}

		### Calculate uncertainty of sample classification into groups
		my $mean_uncertainty   = mean(@{$uncertainty[$n]});
		my $median_uncertainty = median(@{$uncertainty[$n]});
		my $max_uncertainty    = max(@{$uncertainty[$n]});
		my $stddev_uncertainty = stddev(@{$uncertainty[$n]});
		#print $n ."\t". $mean[$n][0] ."\t". $mean[$n][1] ."\t". "$sum_rpkm\t$mean_uncertainty\t$median_uncertainty\t$max_uncertainty\t$stddev_uncertainty\n";


		### Check if bait data is meaningful
		my $use_bait = 1;

		if( $sum_rpkm < 0.01 )              { $use_bait = 0; }   # coverage check
		elsif( $mean_uncertainty > 0.2 )    { $use_bait = 0; }   # uncertainty check
		elsif( $median_uncertainty > 0.2 )  { $use_bait = 0; }   # uncertainty check
		elsif( $proportions[$n][0] < 0.05 ) { $use_bait = 0; }   # minor allele frequency
		elsif( $proportions[$n][0] > 0.95 ) { $use_bait = 0; }   # minor allele frequency
		elsif( ($mean[$n][0] != 0) && ($mean[$n][1]/$mean[$n][0] < $foldchange) ) { $use_bait = 0; }
		

		# All criteria OK
		if( $use_bait == 1 ) {


			# Calculate pattern similarity
			# TODO run through all samples at the last bait ($lastid) and the new bait ($n)

			# Calculate current fold-change
			my $current_fc = 1;
			if($mean[$n][0] == 0) {
				$current_fc = 2;
			}
			else {
				$current_fc = $mean[$n][1]/$mean[$n][0];
			}

			# Start new cluster
			if( $lastid < $n - 10 ) {
				$cluster{$n} = $n;
				$cluster_fc{$n} = $current_fc;
				$cluster_uc{$n} = $mean_uncertainty;
				$cluster_start_id = $n;
				$lastid = $n;
			}
			# Elongagte old cluster
			else {
				$cluster{$cluster_start_id} .= ",$n";
				$cluster_fc{$cluster_start_id} .= ",$current_fc";
				$cluster_uc{$cluster_start_id} .= ",$mean_uncertainty";
				$lastid = $n;
			}

			#print "$n\t$mean[$n][0]\t$mean[$n][1]\n";
		}
	}
}

foreach my $id (sort {$a<=>$b} keys %cluster) {

	my @a = split(",", $cluster{$id});
	my $used_baits = $#a + 1;
	my $baits_in_region = $a[$#a] - $a[0] + 1;
	my $bait_ratio = $used_baits / $baits_in_region;

	if( ($used_baits >= 5)  && ($bait_ratio >= 0.25) ) {

		### Length and used baits per megabase
		my $len = $info[$a[$#a]][1] - $info[$a[0]][0] + 1;
		my $bpm = $used_baits / $len * 1000000;

		### mean and median of fold-change
		my @b = split(",", $cluster_fc{$id});
		my $mean_fc = mean(@b);
		my $median_fc = median(@b);

		### mean uncertainty
		my @c = split(",", $cluster_uc{$id});
		my $mean_uc = mean(@c);

		### Print chr, start, end, length, consecutive baits,  baits in region, baits %, baits/mb, mean fc, median fc
		# TODO: avg-uncertainty, pattern-similarity
		# TODO: printf("%.3f", 3.1415926535);
		print "$id\t$chr\t". $info[$a[0]][0] ."\t". $info[$a[$#a]][1] ."\t$len\t$used_baits\t$baits_in_region\t";
		printf("%.3f\t", $bait_ratio);
		printf("%.3f\t", $bpm);
		#printf("%.3f\t", $mean_fc); TODO: check, is way to high
		printf("%.3f\t", $median_fc);
		printf("%.3f\n", $mean_uc);
	}
}



### Read command line parameters --------------------------------------------------------------
sub GetCom {

	my @usage = ("$0

Files:
--info              STRING
--classification    STRING   Matrix of allele classification (sample vs. bait)
--uncertainty       STRING   Matrix of uncertainties (sample vs. bait)
--means             STRING   Matrix of mean coverage by allele group
--proportions       STRING   Matrix of allele frequencies by group

Params:
--foldchange     INT      Minimum fold change between groups
--chrom          INT      Chromosome
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "info=s", "classification=s", "uncertainty=s", "means=s", "proportions=s", "foldchange=s", "chrom=s");
				
        die("Please specify info file\n") unless defined($CMD{info});
        die("Please specify classification file\n") unless defined($CMD{classification});
        die("Please specify uncertainty file\n") unless defined($CMD{uncertainty});
	die("Please specify means file\n") unless defined($CMD{means});
	die("Please specify proportions file\n") unless defined($CMD{proportions});
	die("Please specify foldchange\n") unless defined($CMD{foldchange});
	die("Please specify chr\n") unless defined($CMD{chrom});

        $file_info           = $CMD{info};
	$file_classification = $CMD{classification};
	$file_uncertainty    = $CMD{uncertainty};
	$file_means          = $CMD{means};
	$file_proportions    = $CMD{proportions};
	$foldchange          = $CMD{foldchange};
	$chr                 = $CMD{chrom};
}

