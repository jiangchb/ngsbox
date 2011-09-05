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
#  Module: Annotation::Genotyping::NGS_Array_Sensitivity.pl
#  Purpose:
#  In:
#  Out:
#


####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################

use DBI;

my $usage = "\n\n$0 database genotype_table shore_table mpileup_table gatk_table\n\n";
my $database       = shift or die $usage;
my $genotype_table = shift or die $usage;
my $shore_table    = shift or die $usage;
my $mpileup_table  = shift or die $usage;
my $gatk_table     = shift or die $usage;


my $dbh = '';
&connect_to_db();


# Get enriched regions
my %enrich = ();
my $e= "SELECT chr, start, end
	FROM exome_enrichment_oligos
	ORDER BY chr, start";
my $sth_e = $dbh->prepare($e);
$sth_e->execute();
while (my $ref = $sth_e->fetchrow_hashref()) {
	my $chr = $ref->{'chr'};
	my $beg = $ref->{'start'};
	my $end = $ref->{'end'};
	
	for (my $i = $beg; $i <= $end; $i++) {
		$enrich{$chr}{$i} = 1;
	}
}


# Get all sample names
my $q= "SELECT distinct sample
	FROM $genotype_table
	ORDER BY sample";
my $sth = $dbh->prepare($q);
$sth->execute();


# For each sample do
while (my $ref = $sth->fetchrow_hashref()) {

	my $sample = $ref->{sample};
	my %shore_results   = ();
	my %mpileup_results = ();
	my %gatk_results    = ();


	# Get all NGS SNP calls from SHORE
	my $q2="SELECT  chr, position, reference, variant, quality, support, frequency
		FROM	$shore_table
		WHERE 	sample = '$sample' &&
			support >= 3 &&
			variant != '-' &&
			frequency >= 0.3 &&
			repetitiveness <= 1.1
		ORDER BY chr, position";
	my $sth2 = $dbh->prepare($q2);
	$sth2->execute();

	while( my $ref2 = $sth2->fetchrow_hashref() ) {
		my $chr = $ref2->{'chr'};
		my $pos = $ref2->{position};

		$shore_results{"$chr-$pos"} = $ref2->{reference} ."\t". $ref2->{variant} ."\t". $ref2->{quality} ."\t". $ref2->{support} ."\t". sprintf("%.3f", $ref2->{frequency});
	}


	# Get all NGS SNP calls from mpileup
	my $q3="SELECT  chr, position, reference, variant, Qs, support, fA, fT, fC, fG
		FROM    $mpileup_table
		WHERE   sample = '$sample' &&
			support >= 3
		ORDER BY chr, position";
	my $sth3 = $dbh->prepare($q3);
	$sth3->execute();

	while (my $ref3 = $sth3->fetchrow_hashref()) {
		my $chr = $ref3->{'chr'};
		my $pos = $ref3->{position};

		my $freq = $ref3->{fA};
		if($ref3->{variant} eq "T") { $freq = $ref3->{fT}; }
		if($ref3->{variant} eq "C") { $freq = $ref3->{fC}; }
		if($ref3->{variant} eq "G") { $freq = $ref3->{fG}; }
		my $support = $ref3->{support} * $freq;

		if($support >= 3 && $freq >= 0.3) {
			$mpileup_results{"$chr-$pos"} = $ref3->{reference} ."\t". $ref3->{variant} ."\t". $ref3->{Qs} ."\t". $support ."\t". sprintf("%.3f", $freq);
		}
	}

	# Get all NGS SNP calls from GATK
	my $q4="SELECT  chr, position, reference, variant, quality, support, frequency
		FROM    $gatk_table
		WHERE   sample = '$sample' &&
			support >= 5
		ORDER BY chr, position";
	my $sth4 = $dbh->prepare($q4);
	$sth4->execute();

	while (my $ref4 = $sth4->fetchrow_hashref()) {
		my $chr = $ref4->{'chr'};
		my $pos = $ref4->{position};

		$gatk_results{"$chr-$pos"} = $ref4->{reference} ."\t". $ref4->{variant} ."\t". $ref4->{quality} ."\t". $ref4->{support} ."\t". sprintf("%.3f", $ref4->{frequency});
	}


	# Get all genotyped positions and print NGS call results
	my $q5="SELECT  *
		FROM    $genotype_table
		WHERE   sample = '$sample' &&
			gtype != 'NC'
		ORDER BY chr, hg19";
	my $sth5 = $dbh->prepare($q5);
	$sth5->execute();

	my %stats = ();

	while (my $ref5 = $sth5->fetchrow_hashref()) {
		my $chr = $ref5->{'chr'};
		my $pos = $ref5->{hg19};

		# Check if part of enriched regions
		if(exists $enrich{$chr}{$pos}) {

			$stats{total_snps}++;
			if($ref5->{gtype} eq "BB") { $stats{BB}++; }
			elsif($ref5->{gtype} eq "AB") { $stats{AB}++; }
			elsif($ref5->{gtype} eq "AA") { $stats{AA}++; }     


			if( exists $shore_results{"$chr-$pos"} || exists $mpileup_results{"$chr-$pos"} ) { #|| exists $gatk_results{"$chr-$pos"} ) {

				### Print Genotyping array result
				print "$sample\t" . $ref5->{name} . "\t$chr\t$pos\t" . $ref5->{gtype} ."\t". sprintf("%.3f", $ref5->{B_allele_freq});

		
				### Print SHORE result
				if( exists $shore_results{"$chr-$pos"} ) {
					my @a = split("\t", $shore_results{"$chr-$pos"});

					print "\t". $shore_results{"$chr-$pos"};
				
					if($ref5->{gtype} eq "AA") {
						$stats{shoreAA}++;
						if(exists $mpileup_results{"$chr-$pos"} ) {
							$stats{shore_mpileupAA}++;
						}
					}
					elsif($ref5->{gtype} eq "AB") {
						$stats{shoreAB}++;
						if(exists $mpileup_results{"$chr-$pos"} ) {
							$stats{shore_mpileupAB}++;
						}
						if($a[4] >= 0.4 && $a[4] <= 0.6) {
							$stats{shoreAB_good}++;
						}
						else {
							$stats{shoreAB_bad}++;
							}
					}
					elsif($ref5->{gtype} eq "BB") {
						$stats{shoreBB}++;
						if(exists $mpileup_results{"$chr-$pos"} ) {
							$stats{shore_mpileupBB}++;
						}
						if($a[4] >= 0.8) {
							$stats{shoreBB_good}++;
						}
						else {
							$stats{shoreBB_bad}++;
						}
					}
				}
				else {
					print "\tNA\tNA\tNA\tNA";
				}


				### Print mpileup result
				if( exists $mpileup_results{"$chr-$pos"} ) {
					my @a = split("\t", $mpileup_results{"$chr-$pos"});

					print "\t". $mpileup_results{"$chr-$pos"};

					if($ref5->{gtype} eq "AA") {
						$stats{mpileupAA}++;
					}
					elsif($ref5->{gtype} eq "AB") {
						$stats{mpileupAB}++;
						if($a[4] >= 0.4 && $a[4] <= 0.6) {
							$stats{mpileupAB_good}++;
						}
						else {
							$stats{mpileupAB_bad}++;
						}
					}
					elsif($ref5->{gtype} eq "BB") {
						$stats{mpileupBB}++;
						if($a[4] >= 0.8) {
							$stats{mpileupBB_good}++;
						}
						else {
							$stats{mpileupBB_bad}++;
						}
					}
				}
				else {
					print "\tNA\tNA\tNA\tNA";
				}
		

				# Print GATK result
				#if( exists $gatk_results{"$chr-$pos"} ) {
				#	print "\t". $gatk_results{"$chr-$pos"};
				#}
				#else {
				#	print "\tNA\tNA\tNA\tNA";
				#}

				print "\n";
			}
		}
	}

	print STDERR "Sample: $sample\n" .
		"Total genotypes: " . $stats{total_snps} .
		"\nAA: " . $stats{AA} . 
		"\nAB: " . $stats{AB} .
		"\nBB: " . $stats{BB} .

		"\n\nSHORE AA: " . $stats{shoreAA} .
		"\nSHORE AB: " . $stats{shoreAB} .
		"\n    good: " . $stats{shoreAB_good} .
		"\n     bad: " . $stats{shoreAB_bad} .
		"\nSHORE BB: " . $stats{shoreBB} .
		"\n    good: " . $stats{shoreBB_good} .
		"\n     bad: " . $stats{shoreBB_bad} .

		"\n\nmpileup AA: " . $stats{mpileupAA} .
		"\nmpileup AB: " . $stats{mpileupAB} .
		"\n      good: " . $stats{mpileupAB_good} .
		"\n       bad: " . $stats{mpileupAB_bad} .
		"\nmpileup BB: " . $stats{mpileupBB} .
		"\n      good: " . $stats{mpileupBB_good} .
		"\n       bad: " . $stats{mpileupBB_bad} .

		"\n\nshore+mpileup AA: " . $stats{shore_mpileupAA} .
		"\nshore+mpileup AB: " . $stats{shore_mpileupAB} .
		"\nshore+mpileup BB: " . $stats{shore_mpileupBB} .
		"\n\n";
}

#####################################################
### Connects to a database and returns databaseHandle
#####################################################
sub connect_to_db
{
	my $databaseName = "$database";
        my $driver = "mysql";
	my $host = "pou.crg.es";
	my $username = "gdvisitor";
	my $password = "";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to genopheno-db. Connect error: $DBI::errstr";
}
