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
#  Module: Analysis::SV::Simulate_SV::check_simulated_insertions.pl
#  Purpose:
#  In:
#  Out:
#




my $usage = "\n$0 simulatedInsertions SVprediction extension\n";

my $insertion = shift or die $usage;
my $prediction = shift or die $usage;
my $extension = shift or die $usage;

#########################################################
# Read in inserted/predicted positions
my %INSPOS_SIM = ();
my %INSPOS_PRED = ();

open FILE, $insertion or die $usage;
while (<FILE>) {
	my @a = split " ";	
	$INSPOS_SIM{$a[0]."#".$a[1]} = $a[3];
}
close FILE;

open FILE, $prediction or die $usage;
while (<FILE>) {
        my @a = split " ";
	for (my $i = $a[4]; $i <= $a[5]; $i++) {
		$INSPOS_PRED{$a[3]."#".$i} = $a[6];
	}
}
close FILE;
#
##########################################################

my $FP = 0; my $FN = 0; my $TP = 0;
my $FP_1 = 0; my $FN_1 = 0; my $TP_1 = 0;
my $FP_11 = 0; my $FN_11 = 0; my $TP_11 = 0; 
my $FP_21 = 0; my $FN_21 = 0; my $TP_21 = 0; 
my $FP_51 = 0; my $FN_51 = 0; my $TP_51 = 0; 
my $FP_101 = 0; my $FN_101 = 0; my $TP_101 = 0; 
my $FP_201 = 0; my $FN_201 = 0; my $TP_201 = 0; 

my %FNstring = ();

##########################################################
# Check deletions
# to assess false negative

open FILE, $insertion or die $usage;
while (<FILE>) {
        my @a = split " ";
	my $len = $a[3];
	my $pos_found = 0;
	for (my $i = $a[1]-$extension; $i<= $a[2]+$extension; $i++) {
		if (defined($INSPOS_PRED{$a[0]."#".$i})) {
			$pos_found++;
		}
	}
	if ($pos_found == 0) {
		$FN++;
		if ($len >= 0 && $len <= 10) {
			$FN_1++;
			$FNstring{1} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
		}
		if ($len >= 11 && $len <= 20) {
                        $FN_11++;
			$FNstring{11} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
                }
		if ($len >= 21 && $len <= 50) {
                        $FN_21++;
			$FNstring{21} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
                }
		if ($len >= 51 && $len <= 100) {
                        $FN_51++;
			$FNstring{51} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
                }
		if ($len >= 101 && $len <= 200) {
                        $FN_101++;
			$FNstring{101} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
                }
		if ($len >= 201 && $len <= 1000) {
                        $FN_201++;
			$FNstring{201} .= $a[0]."\t".$a[1]."\t".$a[2]."\n";
                }
	}
}
close FILE;

#print "######################################\n";
#print "False negative predictions:\n";
#print "#############\n";
#print "Length 11 - 20 bp\n", $FNstring{11};
#print "#############\n";
#print "Length 21 - 50 bp\n", $FNstring{21};
#print "#############\n";
#print "Length 51 - 100 bp\n", $FNstring{51};
#print "#############\n";
#print "Length 101 - 200 bp\n", $FNstring{101};
#print "######################################\n";

#############################################################
# check predictions
# to assess true and false positives

my %FPstring = ();

open FILE, $prediction or die $usage;
open PRED2, ">$prediction.label" or die $usage;

while (my $line = <FILE>) {
	chomp($line);
        my @a = split " ", $line;
	my $len = $a[6];
	my $pos_found = 0;
	for (my $i = $a[4]-$extension; $i<= $a[5]+$extension; $i++) {
		if (defined($INSPOS_SIM{$a[3]."#".$i})) {
			$pos_found++;
		}
	}
	print PRED2 $line;
	if ($pos_found == 0) {
                $FP++;
		print PRED2 "\tF\n"; 
	}
	else {
		$TP++;
		print PRED2 "\tT\n"; 
	}
	if ($len >= 0 && $len <= 10) {
		if ($pos_found == 0) {
        	        $FP_1++;
			$FPstring{1} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_1++;
		}
        }
	if ($len >= 11 && $len <= 20) {
		if ($pos_found == 0) {
        	        $FP_11++;
			$FPstring{11} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_11++;
	        }
	}
	if ($len >= 21 && $len <= 50) {
		if ($pos_found == 0) {
        	        $FP_21++;
			$FPstring{21} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_21++;
		}
        }
	if ($len >= 51 && $len <= 100) {
		if ($pos_found == 0) {
        	        $FP_51++;
			$FPstring{51} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_51++;
		}
        }
	if ($len >= 101 && $len <= 200) {
		if ($pos_found == 0) {
        	        $FP_101++;
			$FPstring{101} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_101++;
		}
        }
	if ($len >= 201) {
		if ($pos_found == 0) {
        	        $FP_201++;
			$FPstring{201} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_201++;
		}
        }

}


#print "######################################\n";
#print "False positive predictions:\n";
#print "#########################\n";
#print "Length 1 - 10 bp\n", $FPstring{1};
#print "Length 11 - 20 bp\n", $FPstring{11};
#print "Length 21 - 50 bp\n", $FPstring{21};
#print "Length 51 - 100 bp\n", $FPstring{51};
#print "Length 101 - 200 bp\n", $FPstring{101};
#print "Length 201 - 1000 bp\n", $FPstring{201};
#print "######################################\n";



my $ppv_1 = 0;
if ($TP_1+$FP_1 != 0) {
	$ppv_1 = $TP_1/($TP_1+$FP_1);
}
my $ppv_11 = 0;
if ($TP_11+$FP_11 != 0) {
        $ppv_11 = $TP_11/($TP_11+$FP_11);
}
my $ppv_21 = 0;
if ($TP_21+$FP_21 != 0) {
        $ppv_21 = $TP_21/($TP_21+$FP_21);
}
my $ppv_51 = 0;
if ($TP_51+$FP_51 != 0) {
        $ppv_51 = $TP_51/($TP_51+$FP_51);
}
my $ppv_101 = 0;
if ($TP_101+$FP_101 != 0) {
        $ppv_101 = $TP_101/($TP_101+$FP_101);
}
my $ppv_201 = 0;
if ($TP_201+$FP_201 != 0) {
        $ppv_201 = $TP_201/($TP_201+$FP_201);
}
my $recall_1 = 0;
if ($TP_1+$FN_1 != 0) {
        $recall_1 = $TP_1/($TP_1+$FN_1);
}
my $recall_11 = 0;
if ($TP_11+$FN_11 != 0) {
        $recall_11 = $TP_11/($TP_11+$FN_11);
}
my $recall_21 = 0;
if ($TP_21+$FN_21 != 0) {
        $recall_21 = $TP_21/($TP_21+$FN_21);
}
my $recall_51 = 0;
if ($TP_51+$FN_51 != 0) {
        $recall_51 = $TP_51/($TP_51+$FN_51);
}
my $recall_101 = 0;
if ($TP_101+$FN_101 != 0) {
        $recall_101 = $TP_101/($TP_101+$FN_101);
}
my $recall_201 = 0;
if ($TP_201+$FN_201 != 0) {
        $recall_201 = $TP_201/($TP_201+$FN_201);
}





print "\n";
print "Length\t\tTP\tFP\tFN\tPPV\tRecall\n";
print "All:\t\t", $TP, "\t", $FP, "\t", $FN, "\t", $TP/($TP+$FP), "\t", $TP/($TP+$FN),"\n";
print "\n";
print "1-10\t\t", $TP_1, "\t", $FP_1, "\t", $FN_1, "\t", $ppv_1, "\t", $recall_1,"\n";
print "11-20\t\t", $TP_11, "\t", $FP_11, "\t", $FN_11, "\t", $ppv_11, "\t", $recall_11,"\n";
print "21-50\t\t", $TP_21, "\t", $FP_21, "\t", $FN_21, "\t", $ppv_21, "\t", $recall_21,"\n";
print "51-100\t\t", $TP_51, "\t", $FP_51, "\t", $FN_51, "\t", $ppv_51, "\t", $recall_51,"\n";
print "101-200\t\t", $TP_101, "\t", $FP_101, "\t", $FN_101, "\t", $ppv_101, "\t", $recall_101,"\n";
print "201-1000\t", $TP_201, "\t", $FP_201, "\t", $FN_201, "\t", $ppv_201, "\t", $recall_201, "\n";




