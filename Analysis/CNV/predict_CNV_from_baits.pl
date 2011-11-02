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


use Getopt::Long;
use FindBin;


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
$i = 0;
my @uncertainty = ();
open UNCERT, $file_uncertainty or die "Cannot open uncertainty file\n";
<UNCERT>;
while( <UNCERT> ) {
	chomp;
	my @a = split("\t", $_);
	$uncertainty[$i] = \@a;
	$i++;
}
close UNCERT or die;

#-> Sum up uncertainty per bait. avg must be < 0.001


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
my $cluster_start_id = -10;
my $lastid = -10;
my $conseq = 0; # is this useful? maybe set to 0 if patterns do not overlap or if only one class was predicted

for(my $n = 0; $n <= $#mean; $n++) {

	# Criteria 1: at least two groups
	if( ($mean[$n][1] ne "NA") && ($mean[$n][1] >= 0.01 ) ) {

		# Criteria 2: Minor allele frequency
		if( ($proportions[$n][0] >= 0.03) && ($proportions[$n][1] >= 0.03) ) {

			# Criteria 3: fold change
			if( ($mean[$n][0] == 0) || ($mean[$n][1]/$mean[$n][0] >= $foldchange) ) {

				# Calculate pattern similarity
				# TODO run through all samples at the last bait ($lastid) and the new bait ($n)

				# Start new cluster
				if( $lastid < $n - 10 ) {
					$cluster{$n} = $n;
					$cluster_start_id = $n;
					$lastid = $n;
				}
				# Elongagte old cluster
				else {
					$cluster{$cluster_start_id} .= ",$n";
					$lastid = $n;
				}

				#print "$n\t$mean[$n][0]\t$mean[$n][1]\n";
			}
		}
	}
}

foreach my $id (sort {$a<=>$b} keys %cluster) {
	my @a = split(",", $cluster{$id});
	my $count = $#a + 1;
	my $bait_ratio = $#a / ($a[$#a] - $a[0] + 1);

	if( ($count >= 4)  && ($bait_ratio > 0.25) ) {
		#print "COUNT: " . $count . "\n";
		#print $cluster{$id} . "\n";
		my $len = $info[$a[$#a]][1] - $info[$a[0]][0] + 1;

		### Print chr, start, end, length, consecutive baits, catched baits %, avg-uncertainty, pattern-similarity, avg-foldchange
		print $chr ."\t". $info[$a[0]][0] ."\t". $info[$a[$#a]][1] ."\t$len\t$count\t$bait_ratio\n";
	}
}

### Min/Max calculation
sub min {
	my ($a, $b) = @_;
	return $a if $a < $b;
	return $b;
}
sub max {
	my ($a, $b) = @_;
	return $a if $a > $b;
	return $b;
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

