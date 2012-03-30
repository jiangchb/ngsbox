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
#  Module: Analysis::Assembly::ImproveAssembly::improve_assembly_by_realignment.pl
#  Purpose:
#  In:
#  Out:
#



### User params
my $usage = "\n$0 trim_len min_length min_read_per_ctg min_bp_per_mm max_obs_exp_cov contig_file config_file unseq_file snp_file deletion_file insertion_file\n\n";

my $trim_len         = shift or die $usage;	# Number of bases that are trimmed from the 5' and 3' end of a contig by default
my $min_length       = shift or die $usage;	# Filter threashold: Minimum length of a contig
my $min_read_per_ctg = shift or die $usage;	# Filter threashold: Minimum number of reads aligning to a contig
my $min_bp_per_mm    = shift or die $usage;	# Filter threashold: mismatches per x based (e.g. for 1 error per 10kb define min_bp_per_mm=10000)
my $max_obs_exp_cov  = shift or die $usage;	# Filter threashold: Maximum ratio of observed to expected coverage

my $contig_file    = shift or die $usage;	# fasta file with assembled contigs
my $config_file    = shift or die $usage;	# config file of 'shore consensus'
my $unseq_file     = shift or die $usage;	# unsequenced file of 'shore consensus'
my $snp_file       = shift or die $usage;	# SNP file of 'shore consensus'
my $deletion_file  = shift or die $usage;	# Deletion file of 'shore consensus'
my $insertion_file = shift or die $usage;	# Insertion file of 'shore consensus'


### Init containers
my %ctg_seq              = ();
my %ctg_len              = ();
my %ctg_readcount        = ();
my %ctg_unique_readcount = ();
my %ctg_depth            = ();
my %ctg_unique_depth     = ();

my %splitter = ();
my %unseq_N  = ();
my %snp      = ();

my $id = -1;


### Parse contigs
open CONTIG, $contig_file or die "Cannot open $contig_file\n";
while(<CONTIG>) {
	chomp($_);

	if (substr($_, 0, 1) eq ">") {
		$id = substr($_, 1);
		$id =~ s/scaffold_//g;
	}
	else {
		$ctg_seq{$id} .= $_;
	}
}
close CONTIG;


### Set contig length
my $total_contig_len = 0;
foreach $id (keys %ctg_seq) {
	$ctg_len{$id} = length($ctg_seq{$id});
	$total_contig_len += $ctg_len{$id};
}


### Read re-alignment stats
my $genome_unique_depth = 1;
open CONFIG, $config_file or die "Cannot open $config_file\n";
while(<CONFIG>) {

	chomp($_);

	if(substr($_, 0, 12) eq "DEPTH GENOME") {
		my ($junk, $value) = split(": ", $_);
		my($unique_bases_sequenced, $unique_positions, $sequencing_depth) = split(/ .{1} /, $value);
		$genome_unique_depth = $sequencing_depth;
	}

	elsif(substr($_, 0, 9) eq "DEPTH CHR") {
		
		# get contig id
		my ($id, $value) = split(": ", $_);
		$id =~ s/DEPTH CHR-//g;
		$id--;
		
		# split values
		my($unique_bases_sequenced, $unique_positions, $unique_depth) = split(/ .{1} /, $value);
		$ctg_unique_depth{$id} = $unique_depth;
	}

	elsif(substr($_, 0, 16) eq "UNIQUE READS CHR") {
		
		# get contig id
		my ($id, $readcount) = split(": ", $_);
		$id =~ s/UNIQUE READS CHR-//g;
		$id--;
		$ctg_unique_readcount{$id} = $readcount;
	}

	elsif(substr($_, 0, 15) eq "TOTAL READS CHR") {
		my ($id, $readcount) = split(": ", $_);
		$id =~ s/TOTAL READS CHR-//g;
		$id--;
		$ctg_readcount{$id} = $readcount;
	}
}


### Set avg depth (including rep reads
my $total_yield  = 0;
foreach $id (keys %ctg_seq) {
	$ctg_depth{$id} = $ctg_readcount{$id} * 100 / $ctg_len{$id};
	$total_yield += $ctg_readcount{$id} * 100;
}
my $genome_depth = $total_yield / $total_contig_len;


### Parse unsequenced regions
open UNSEQ, $unseq_file or die "Cannot open $unseq_file\n";
while( <UNSEQ> ) {

	my @a = split("\t", $_);
	my $id = $a[1];
	$id =~ s/scaffold_//g;

	if(! exists $unseq_N{$id} ) {
		my @tmp1 = ();
		my @tmp2 = ();
		$unseq_N{$id} = \@tmp1;
		$splitter{$id} = \@tmp2;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		my $unseq_seq = substr($ctg_seq{$id}, $a[2] - 1, $a[3] - $a[2] + 1);

		if($unseq_seq =~ /N/) {
			push( @{$unseq_N{$id}}, "$a[2]-$a[3]" );
		}
		else {
			push( @{$splitter{$id}}, "$a[2]-$a[3]" );
		}
	}
}


### Parse SNPs
open SNP, $snp_file or die "Cannot open $snp_file\n";
while( <SNP> ) {
	my @a = split("\t", $_);
	my $id = $a[1];
	$id =~ s/scaffold_//g;

	if(! exists $snp{$id} ) {
		my %tmp = ();
		$snp{$id} = \%tmp;
	}

	# Do not count in trimming range
	if( ($a[2] > $trim_len) && ($a[2] < ($ctg_len{$id} - $trim_len)) ) {

		# check IUPAC code of contig
		if( ($a[3] eq "A" || $a[3] eq "C" || $a[3] eq "G" || $a[3] eq "T") && ($a[4] ne "-") ) {
		
			# Select high quality SNPs 
			if($a[5] >= 30 && $a[7] >= 0.7) {
				$snp{$id}{$a[2]} = 1;
			}
		}
	}
}


### Parse small deletions
open DELETION, $deletion_file or die "Cannot open $deletion_file\n";
while( <DELETION> ) {
        my @a = split("\t", $_);
	my $id = $a[1];
	$id =~ s/scaffold_//g;

	if(! exists $splitter{$id} ) {
		my @tmp = ();
		$splitter{$id} = \@tmp;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		push( @{$splitter{$id}}, "$a[2]-$a[3]");
	}
}


### Parse small insertions
open INSERTION, $insertion_file or die "Cannot open $insertion_file\n";
while( <INSERTION> ) {
        my @a = split("\t", $_);
	my $id = $a[1];
	$id =~ s/scaffold_//g;

	if(! exists $splitter{$id} ) {
		my @tmp = ();
		$splitter{$id} = \@tmp;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		push( @{$splitter{$id}}, "$a[2]-$a[3]");
	}
}


# Mask sequence around gaps
foreach $id ( keys %splitter ) {
	my %split_locus = ();

	foreach my $locus ( @{$splitter{$id}} ) {
		my ($s, $e) = split("-", $locus);
		$s -= 20; $e += 20;

		for(my $i = $s; $i <= $e; $i++) {
			substr($ctg_seq{$id}, $i-1, 1, "N");
			push( @{$unseq_N{$id}}, "$s-$e" );
			delete $snp{$id}{$i};
		}
	}
}


### Validate
foreach $id (sort {$a<=>$b} keys %ctg_seq) {
	
	my $remapping_mm =  scalar( keys %{$snp{$id}} );

	# Remove low read count contigs
	if( $ctg_readcount{$id} < $min_read_per_ctg) {
		delete $ctg_seq{$id};
		print STDERR "BAD COUNT: $id\t" . $ctg_len{$id} . "\t" . $ctg_readcount{$id} . "\n";
	}

	# Remove bad coverage contigs
	elsif(	($ctg_unique_depth{$id} > 0) &&
		(
		  (($ctg_unique_depth{$id} / $genome_unique_depth) >= $max_obs_exp_cov) ||
		  #(($genome_unique_depth / $ctg_unique_depth{$id}) >= $max_obs_exp_cov)
		  ($ctg_unique_depth{$id} < 5)
		)
	){
		delete $ctg_seq{$id};
		print STDERR "BAD COVERAGE: $id\t" . $ctg_len{$id} . "\t" . $ctg_depth{$id} . "\n";
	}

	# Remove erroneous contigs
	elsif( ($remapping_mm != 0) && (($ctg_len{$id} / $remapping_mm) < $min_bp_per_mm) ) {
		delete $ctg_seq{$id};
		print STDERR "MISMATCHES: $id\t" . $ctg_len{$id} . "\t" . $ctg_len{$id} / $remapping_mm . "\n";
	}

	# Remove short contigs
	elsif($ctg_len{$id} < $min_length) {
		delete $ctg_seq{$id};
		print STDERR "BAD LENGTH: $id\t" . $ctg_len{$id} . "\n";
	}
}


### Print validated contigs
my $total_genome_coverage = 0;
foreach $id (sort {$a<=>$b} keys %ctg_seq) {
	
	$ctg_len{$id} -= (2 * $trim_len);
	$total_genome_coverage += $ctg_len{$id};
	my $trimmed_contig = substr($ctg_seq{$id}, $trim_len, $ctg_len{$id});

	print ">$id | " . $ctg_len{$id} . "\n" . $trimmed_contig . "\n";
}

print STDERR "Average contig coverage: $genome_depth\n";
print STDERR "Total genome covered: $total_genome_coverage\n";


exit(0);
