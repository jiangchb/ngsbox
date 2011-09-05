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
#  Module: Simulation::Simulate_NGS::random_insertions_by_deletions.pl
#  Purpose:
#  In:
#  Out:
#


##############################################
# Simulates insertions based on a empirical 
# distribution of deletion lengths from the
# same strain.
##############################################

our $VERSION = '1.0';


use Getopt::Long;

# User parameters
my $deletion_file     = "";
my $genome1_file      = "";
my $genome2_file      = "";
my $transposon_file   = "";
my $annotation_file   = "";
my $out_file          = "";
my $transposon_copies = 0;
my $chr_len           = 30432563;

my %CMD;
GetCom();


# Hashes and arrays
my @nucleotides = ('A', 'T', 'C', 'G');
my %transposons = ();
my %del = ();
my %ins = ();


### Read sequence from genome 1
open( GENOME1, $genome1_file ) or die "Cannot open input file $genome1_file\n";
my $genome_seq_1 = "";
while( <GENOME1> ) {
	chomp;
	if($_ !~ />/) {
		$_ =~ s/N//g;
		$genome_seq_1 .= uc($_);
	}
}


### Read sequence from genome 2
open GENOME2, $genome2_file or die "Cannot open input file $genome2_file\n";
my $genome_seq_2 = "";
while( <GENOME2> ) {
	chomp;
	if($_ !~ />/) {
		$_ =~ s/N//g;
		$genome_seq_2 .= uc($_);
	}
}


### Read fasta file of transposable elements (third source for insertion sequence
open TE, $transposon_file or die "Cannot open input file $transposon_file\n";
my $te_seq = "";
my $te_count = -1;
while( <TE> ) {
	chomp;
	if($_ =~ />/) {
		if($te_count != -1) {
			$transposons{$te_count} = $te_seq;
			$te_seq = "";
		}
		$te_count++;
	}
	else {
		$te_seq .= uc($_);
	}
}


### Read deletion file
open DELETION, $deletion_file or die "Cannot open file $deletion_file\n";
while( <DELETION> ) {
	chomp;
	my ($beg, $end) = split("\t", $_);
	$del{$beg} = $end - $beg + 1;
}


### Create random set of insertion based on the deletion distribution
foreach my $pos ( sort {$a<=>$b}  keys %del ) {
	my $len = $del{$pos};
	
	# Randomize position
	my $select_genome = int(rand(2));
	my $ins_pos = int(rand($chr_len)) + 1;
	while( exists $ins{$ins_pos}) {
		$ins_pos = int(rand($chr_len)) + 1;
	}

	if( ! exists $ins{$ins_pos} ) {

		# Short insertions have random seq
		if($len <= 3) {
			my $ins_event = int(rand(4));
			$ins{$ins_pos} = $nucleotides[$ins_event];

			if($len > 1) {
				my $ins_event = int(rand(4));
				$ins{$ins_pos} .= $nucleotides[$ins_event];
			}
				
			if($len > 2) {
				my $ins_event = int(rand(4));
				$ins{$ins_pos} .= $nucleotides[$ins_event];
			}
		}
		else {
			# Insert long seq from genome 1
			if($select_genome == 0) {
				my $random_pos = int(rand(length($genome_seq_1)));
				$ins{$ins_pos} = substr($genome_seq_1, $random_pos, $len);
			}

			# Insert long seq from genome 2
			elsif($select_genome == 1) {
				my $random_pos = int(rand(length($genome_seq_2)));
				$ins{$ins_pos} = substr($genome_seq_2, $random_pos, $len);
			}
		}
	}
}


### Insert 10 random transposons
for( my $i = 0; $i < $transposon_copies; $i++ ) {
	my $transposon_num = (int(rand($te_count)));
	my $ins_pos = int(rand($chr_len)) + 1;
	while( exists $ins{$ins_pos}) {
		$ins_pos = int(rand($chr_len)) + 1;
	}

	$ins{$ins_pos} = $transposons{$transposon_num};
}


### Write insertions to file
open OUT, ">$out_file" or die "Cannot open output file $out_file\n";
foreach my $ins_pos (sort {$a<=>$b} keys %ins) {
	print OUT "$ins_pos\t" . $ins{$ins_pos} . "\n";
}


####################################################
# Randomize new base
sub randomize_base {
	my $base = shift;
	my @nucleotides = ('A', 'T', 'C', 'G');
	my $mutated_base = $base;
	while($mutated_base eq $base) {
		my $snp_event = int(rand(4));
		$mutated_base = $nucleotides[$snp_event];
	}
	return($mutated_base);
}


### Read command line parameters
sub GetCom {

	my @usage = ("$0

--delfile      STRING    Deletions (used to define count and length distribution of simulated insertions)
--genome1      STRING    Chromosome from different species in fasta format
--genome2      STRING    Chromosome from different species in fasta format
--transposon   STRING    Transposons in fasta format
--outfile      STRING    Outfile for simulated insertions

Optional:
--annotation   STRING    Gene annotation in GFF format            (default: none)
--chrlength    INT       Length of reference chromosome           (default: 30,432,563)
--TEcopies     INT       Number of simulated transposon copies    (default: 0)

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "delfile=s","genome1=s", "genome2=s", "transposon=s", "outfile=s", "annotation=s", "chrlength=s", "TEcopies=s");

	die("Please specify deletion infile\n") unless defined($CMD{delfile});
	die("Please specify genome1 infile\n") unless defined($CMD{genome1});
	die("Please specify genome2 infile\n") unless defined($CMD{genome2});
	die("Please specify transposon infile\n") unless defined($CMD{transposon});
	die("Please specify outfile\n") unless defined($CMD{outfile});

	$deletion_file   = $CMD{delfile};
	$genome1_file    = $CMD{genome1};
	$genome2_file    = $CMD{genome2};
	$transposon_file = $CMD{transposon};
	$out_file        = $CMD{outfile};

	if (defined($CMD{annotation})) { $annotation_file = $CMD{annotation}; }
	if (defined($CMD{chrlength}))  { $chr_len = $CMD{chrlength}; }
	if (defined($CMD{TEcopies}))  { $transposon_copies = $CMD{TEcopies}; }
}
