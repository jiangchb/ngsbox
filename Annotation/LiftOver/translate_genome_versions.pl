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
#  Module: Annotation::LiftOver::translate_genome_versions.pl
#  Purpose:
#  In:
#  Out:
#



# translates positions of different reference genome versions.
# optimized for Thermotoga (bacteria)

my $usage = "$0 snpfile deletion_file insertion_file long_insert fasta_file\n";
my $snpfile = shift or die $usage;
my $delfile = shift or die $usage;
my $insfile = shift or die $usage;
my $lifile  = shift or die $usage;
my $fafile  = shift or die $usage;


### Load SNPs
open SNP, $snpfile or die "Cannot open $snpfile\n";
my %SNPS = ();
while( <SNP> ) {
	chomp;
	my @a = split " ", $_;
	$SNPS{$a[2]} = $a[4];
}
close SNP;


### Load Deletions
open DEL, $delfile or die "Cannot open $delfile\n";
my %DELS = ();
while( <DEL> ) {
	chomp;
	my @a = split " ", $_;
	$DELS{$a[2]} = $a[5];
}
close DEL;


### Load Insertions
open INS, $insfile or die "Cannot open $insfile\n";
my %INSS = ();
while( <INS> ) {
	chomp;
	my @a = split " ", $_;
	$INSS{$a[3]} = $a[5];
}


### Load Long Insertion
open LI, $lifile or die "Cannot open $lifile\n";
my %LIS = ();
while( <LI> ) {
	chop;
	my @a = split " ", $_;
	$LIS{$a[2]} = $a[3];
}


### Parse fasta
open FASTA, $fafile or die "Cannot open $fafile\n";
my $seq = "";
while ( <FASTA> ) {
	chomp;
	if (substr($_, 0, 1) ne ">") {
		$seq .= $_;
	}
}


### Translate positions

my $pos_A = 1;
my $pos_B = 1;
my $pos_C = 1;


for (my $i = 0; $i < length($seq); $i++) {

	my $nuc = substr($seq, $i, 1);

	if( defined($LIS{$pos_B} ) ) {

		for(my $j = 0; $j < 8869; $j++) {
			#print STDERR "NA\t-\t$pos_B\t$nuc\t$pos_C\t$nuc\tL\n";
			print "NA\t-\t$pos_B\t$nuc\t$pos_C\t$nuc\tL\n";
			$pos_B++;
			$pos_C++;
			$i++;
			$nuc = substr($seq, $i, 1);
		}
		print "$pos_A\t$nuc\t$pos_B\t$nuc\t$pos_C\t$nuc\tR\n";
		$pos_A++;
		$pos_B++;
		$pos_C++;
	}

	elsif( defined($DELS{$pos_B}) ) {
		print "$pos_A\t$nuc\t$pos_B\t$nuc\tNA\t-\tD\n";
		$pos_A++;
		$pos_B++;
	}

	elsif( defined($INSS{$pos_B}) ) {
		for(my $j = 0; $j < length($INSS{$pos_B}); $j++) {
			print"NA\t-\tNA\t-\t$pos_C\t". substr($INSS{$pos_B}, $j, 1) . "\tI\n";
			$pos_C++;
		}
		print "$pos_A\t$nuc\t$pos_B\t$nuc\t$pos_C\t$nuc\tR\n";
		$pos_A++;
		$pos_B++;
		$pos_C++;
	}

	elsif( defined($SNPS{$pos_B}) ) {
		print "$pos_A\t$nuc\t$pos_B\t$nuc\t$pos_C\t" . $SNPS{$pos_B} ."\tS\n";
		$pos_A++;
		$pos_B++;
		$pos_C++;
	}

	else {
		print "$pos_A\t$nuc\t$pos_B\t$nuc\t$pos_C\t$nuc\tR\n";
		$pos_A++;
		$pos_B++;
		$pos_C++;
	}
}

