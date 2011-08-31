#! /usr/bin/perl

use strict;

my $usage = "$0 map.list snp.txt reference.txt deletion.txt insertion.txt\n";
my $map = shift or die $usage;
my $snp = shift or die $usage;
my $ref = shift or die $usage;
my $del = shift or die $usage;
my $ins = shift or die $usage;
open MAP, $map or die $usage;
open SNP, $snp or die $usage;
open REF, $ref or die $usage;
open DEL, $del or die $usage;
open INS, $ins or die $usage;

my %INSS = ();
while (my $line = <INS>) {
        my @a = split " ", $line;
        $INSS{$a[1]."#".$a[2]} = $a[4];
}
close INS;
print STDERR "Got insertions\n";

my %DELS = ();
while (my $line = <DEL>) {
        my @a = split " ", $line;
	for (my $i = $a[2]; $i <= $a[3]; $i++) {
	        $DELS{$a[1]."#".$i} = 1;
	}
}
close DEL;
print STDERR "Got deletions\n";

my %REFS = ();
while (my $line = <REF>) {
        my @a = split " ", $line;
	if ($a[2] > 5000000) {
		last;
	}
        $REFS{$a[1]."#".$a[2]} = 1;
}
close REF;
print STDERR "Got refernces\n";

my %SNPS = ();
while (my $line = <SNP>) {
	my @a = split " ", $line;
	$SNPS{$a[1]."#".$a[2]} = 1;
}
close SNP;
print STDERR "Got SNPa\n";

my $SNP_MM = 0;
my $RE_MM = 0;

my $DEL_GAP = 0;
my $RE_GAP = 0;

my $INS_INS = 0;
my $RE_INS = 0;

my $count = 0;
RUN: while (my $line = <MAP>) {
	$count++;
	if ($count % 100000 == 0) {
		#print $count, "\n";
	}
	my @a = split " ", $line;
	my $chr = $a[0];
	my $pos = $a[1];
	if ($pos > 5000000) {
		last RUN;
	}
	if ($a[6] == 1) { 
		for (my $i = 0; $i < length($a[2]); $i++) {
			if (substr($a[2], $i, 1) eq "[") {
				if (substr($a[2], $i+1, 1) ne "-" and substr($a[2], $i+2, 1) ne "-") { # no gap
					if (defined($SNPS{$chr."#".$pos})) {
						$SNP_MM++;
					}
					elsif (defined($REFS{$chr."#".$pos})) {
						$RE_MM++;
					}
				}
				if (substr($a[2], $i+1, 1) ne "-" and substr($a[2], $i+2, 1) eq "-") { # gap in read
					if (defined($DELS{$chr."#".$pos})) {
                                                $DEL_GAP++;
                                        }
                                        elsif (defined($REFS{$chr."#".$pos})) {
                                                $RE_GAP++;
                                        }
				}
				if (substr($a[2], $i+1, 1) eq "-" and substr($a[2], $i+2, 1) ne "-") { # gap in ref
					if (defined($INSS{$chr."#".$pos})) {
                                                $INS_INS++;
                                        }
                                        else {
                                                $RE_INS++;
                                        }
                                }
				#print $chr, "\t", $pos, "\t", $SNP_MM, "\t", $RE_MM, "\n";
				if (!(substr($a[2], $i+1, 1) eq "-")) {
					$pos++;
				}
				$i += 3;
			}
			else {
				$pos++;
			}
		}
	}
}


print "# Mismatches due to SNPs: ", $SNP_MM,"\t", $SNP_MM/($SNP_MM+$RE_MM), "\n";
print "# Mismatches due to sequencing errors: ", $RE_MM, "\t", $RE_MM/($SNP_MM+$RE_MM), "\n";
print "\n";
print "# Gaps due to Deletions: ", $DEL_GAP,"\t", $DEL_GAP/($DEL_GAP+$RE_GAP), "\n";
print "# Mismatches due to sequencing errors: ", $RE_GAP, "\t", $RE_GAP/($DEL_GAP+$RE_GAP), "\n";
print "\n";
print "# Insertions due to Inserions: ", $INS_INS,"\n";
print "# Insertions due to sequencing errors: ", $RE_INS, "\n";



