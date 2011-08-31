#!/usr/bin/perl

use strict;


my $usage = "\n$0 simulatedDeletions SVprediction fraction\n";

my $deletion = shift or die $usage;
my $prediction = shift or die $usage;
my $FOUND = shift or die $usage;

#########################################################
# Read in deleted/predicted positions
my %DELPOS_SIM = ();
my %DELPOS_PRED = ();

print STDERR "###########################\n";

open FILE, $deletion or die $usage;
while (<FILE>) {
	my @a = split " ";	
	for (my $i = $a[1]; $i<= $a[2]; $i++) {
		$DELPOS_SIM{$a[0]."#".$i} = 1;
	}
}
close FILE;


print STDERR "###########################\n";

open FILE, $prediction or die $usage;
while (<FILE>) {
        my @a = split " ";
        for (my $i = $a[4]; $i<= $a[5]; $i++) {
                $DELPOS_PRED{$a[3]."#".$i} = 1;
        }
}
close FILE;

##########################################################

my $FP = 0; my $FN = 0; my $TP = 0;
my $FP_1 = 0; my $FN_1 = 0; my $TP_1 = 0;
my $FP_11 = 0; my $FN_11 = 0; my $TP_11 = 0; 
my $FP_21 = 0; my $FN_21 = 0; my $TP_21 = 0; 
my $FP_51 = 0; my $FN_51 = 0; my $TP_51 = 0; 
my $FP_101 = 0; my $FN_101 = 0; my $TP_101 = 0; 
my $FP_201 = 0; my $FN_201 = 0; my $TP_201 = 0; 

my %FNstring = ();

print STDERR "###########################\n";

##########################################################
# Check deletions
# to assess false negative

open FILE, $deletion or die $usage;
while (<FILE>) {
        my @a = split " ";
	my $len = $a[3];
	my $pos_found = 0;
	for (my $i = $a[1]; $i<= $a[2]; $i++) {
		if (defined($DELPOS_PRED{$a[0]."#".$i})) {
			$pos_found++;
		}
	}
	if ($pos_found == 0 || $pos_found/$len < $FOUND) {
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

print STDERR "###########################\n";

#print "######################################\n";
#print "False negative predictions:\n";
##print $FNstring{1}, $FNstring{11}, $FNstring{21}, $FNstring{51}, $FNstring{101}, $FNstring{201};
#print "#############\n";
#print "Length 51 - 100 bp\n", $FNstring{51};
#print "#############\n";
#print "Length 101 - 200 bp\n", $FNstring{101};
#print "#############\n";
#print "Length 201 - 1000 bp ", $FNstring{201};
#print "######################################\n";

#############################################################
# check predictions
# to assess true and false positives

my %FPstring = ();

open FILE, $prediction or die $usage;

while (<FILE>) {
        my @a = split " ";
	my $len = $a[7];
	my $pos_found = 0;
	for (my $i = $a[4]; $i<= $a[5]; $i++) {
		if (defined($DELPOS_SIM{$a[3]."#".$i})) {
			$pos_found++;
		}
	}
	if ($pos_found == 0 || $pos_found/$len < $FOUND) {
                $FP++;
	}
	else {
		$TP++;
	}
	if ($len >= 0 && $len <= 10) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_1++;
			$FPstring{1} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_1++;
		}
        }
	if ($len >= 11 && $len <= 20) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_11++;
			$FPstring{11} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_11++;
	        }
	}
	if ($len >= 21 && $len <= 50) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_21++;
			$FPstring{21} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_21++;
		}
        }
	if ($len >= 51 && $len <= 100) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_51++;
			$FPstring{51} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_51++;
		}
        }
	if ($len >= 101 && $len <= 200) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_101++;
			$FPstring{101} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_101++;
		}
        }
	if ($len >= 201) {
		if (($pos_found == 0 || $pos_found/$len < $FOUND)) {
        	        $FP_201++;
			$FPstring{201} .= $a[3]."\t".$a[4]."\t".$a[5]."\n";
	        }
        	else {
                	$TP_201++;
		}
        }

}

print STDERR "###########################\n";

#print "######################################\n";
#print "False positive predictions:\n";
#print "#########################";
#print "Length 1 - 10 bp\n", $FPstring{1};
#print "Length 11 - 20 bp\n", $FPstring{11};
#print "Length 21 - 50 bp\n", $FPstring{21};
print "Length 51 - 100 bp\n", $FPstring{51};
#print "Length 101 - 200 bp\n", $FPstring{101};
#print "Length 201 - 1000 bp\n", $FPstring{201};
##print $FPstring{51}, $FPstring{101}, $FPstring{201};
#print "######################################\n";


print "\n";
print "Length\tTP\tFP\tFN\tPPV\tRecall\n";
print "All:\t", $TP, "\t", $FP, "\t", $FN, "\t", $TP/($TP+$FP), "\t", $TP/($TP+$FN),"\n";
print "\n";
print "1-10\t", $TP_1, "\t", $FP_1, "\t", $FN_1, "\t", $TP_1/($TP_1+$FP_1), "\t", $TP_1/($TP_1+$FN_1),"\n";
print "11-20\t", $TP_11, "\t", $FP_11, "\t", $FN_11, "\t", $TP_11/($TP_11+$FP_11), "\t", $TP_11/($TP_11+$FN_11),"\n";
print "21-50\t", $TP_21, "\t", $FP_21, "\t", $FN_21, "\t", $TP_21/($TP_21+$FP_21), "\t", $TP_21/($TP_21+$FN_21),"\n";
print "51-100\t", $TP_51, "\t", $FP_51, "\t", $FN_51, "\t", $TP_51/($TP_51+$FP_51), "\t", $TP_51/($TP_51+$FN_51),"\n";
print "101-200\t", $TP_101, "\t", $FP_101, "\t", $FN_101, "\t", $TP_101/($TP_101+$FP_101), "\t", $TP_101/($TP_101+$FN_101),"\n";
print "201-1000\t", $TP_201, "\t", $FP_201, "\t", $FN_201, "\t", $TP_201/($TP_201+$FP_201), "\t", $TP_201/($TP_201+$FN_201),"\n";



