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
#  Module: Analysis::Assembly::Overlap::overlap_control_db.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

### User params
my $usage = "\n\n$0 sample control number_chromosomes min_pairs max_pvalue variants_min_length variants_max_length stddev database table1 table2\n\n";
my $sample     = shift or die $usage;
my $controls   = shift or die $usage;
my $num_chr    = shift or die $usage;
my $min_pairs  = shift or die $usage;
my $max_pvalue = shift or die $usage;
my $min_length = shift or die $usage;
my $max_length = shift or die $usage;
my $stddev      = shift or die $usage;
my $database   = shift or die $usage;
my $table1     = shift or die $usage;
my $table2     = shift or die $usage;

my $dbh;
&connect_to_db();


my %ov_len_chr = ();
my %total_len_chr = ();
my %detected_chr = ();
my %total_chr = ();

my $ov_len_genome = 0;
my $total_len_genome = 0;
my $detected_genome = 0;
my $total_genome = 0;

$stddev *= 2;

for(my $chr = 1; $chr <= $num_chr; $chr++) {

	### Load control samples from table 1
	my @control_ids = split(",", $controls);
	my %control_variants = ();
	my %control_locus = ();

	foreach my $control ( @control_ids ) {

		my $q ="SELECT 	begin, end
			FROM 	$table1
			WHERE 	sample = '$control'
				&& chromosome = $chr
				# && pairs >= $min_pairs
				# && pvalue <= $max_pvalue
				# && p_l_length between $min_length AND $max_length
			ORDER BY begin
		";
		my $sth = $dbh->prepare($q);
		$sth->execute();

		while(my $ref = $sth->fetchrow_hashref()) {

			$control_locus{$ref->{begin}} = $ref->{end};

			for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
				$control_variants{$i} = 1;
			}
		}
	}


	### Load ranged features of table 2 and compare
	my $q ="SELECT 	chromosome, begin, end, pairs, p_l_length, length, pvalue
		FROM	$table2
		WHERE 	sample = '$sample'
			&& chromosome = $chr
			&& pairs >= $min_pairs
			&& pvalue <= $max_pvalue
			&& p_l_length between $min_length AND $max_length
			&& length between $min_length AND $max_length
			&& (p_l_length / length) between 0.5 AND 1.9
		ORDER by chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		
		# Chromosome and geneome wide overlap calculation
		$total_chr{$chr}++;
		$total_genome++;
		$total_len_chr{$chr} += ($ref->{end} - $ref->{begin} + 1);
		$total_len_genome += ($ref->{end} - $ref->{begin} + 1);


		# SV specific overlap
		my $found = 0;
		my $breakpoint_overlap = 0;
		my $SV_overlap_sample_control = 0;
		my $sample_SV_len = $ref->{end} - $ref->{begin} + 1;

		# Check breakpoint overlap
		foreach my $control_begin (keys %control_variants) {
			if( 	($control_begin >= $ref->{begin} - $stddev) && 
				($control_begin <= $ref->{begin} + $stddev) &&
				($control_variants{$control_begin} >= $ref->{end} - $stddev) &&
				($control_variants{$control_begin} <= $ref->{end} + $stddev)
			) {
				$breakpoint_overlap = 1;
			}
		}

		# Check overlap for SV
		for(my $pos = $ref->{begin}; $pos <= $ref->{end}; $pos++) {
			if(exists $control_variants{$pos}) {
				$SV_overlap_sample_control++;
				$ov_len_chr{$chr}++;
				$ov_len_genome++;
				$found = 1;
			}
		}

		# Check if overlapping at all (by at least one base)
		if( $found == 1) {
			$detected_chr{$chr}++;
			$detected_genome++;
		}

		# Print sample specific SVs
		if( ($SV_overlap_sample_control / $sample_SV_len < 0.5) && ($breakpoint_overlap == 0) ) {
			print "$sample\t" . $ref->{chromosome} . "\t" . $ref->{begin} . "\t" . $ref->{end} . "\t" . 
				$ref->{pairs} . "\t" . $ref->{pvalue} . "\t" . $ref->{p_l_length} . "\t". $ref->{'length'} . "\n";
		}
	}

	# Print chromosome wide overlap statistics
	if( ($total_chr{$chr} > 0) && ($total_len_chr{$chr} > 0) ) {
		print STDERR "Detected: Chr$chr\t" . 
			$detected_chr{$chr} ."\t". $total_chr{$chr} ."\t". 
			$detected_chr{$chr}/$total_chr{$chr} . "\n";
		
		print STDERR "Overlap:  Chr$chr\t" . 
			$ov_len_chr{$chr} ."\t". $total_len_chr{$chr} ."\t". 
			$ov_len_chr{$chr}/$total_len_chr{$chr} . "\n\n";
	}
}

# Print Genome wide overlap statistics
if( ($total_genome > 0) && ($total_len_genome > 0) ) {
	print STDERR "Detected: Genome\t$detected_genome\t$total_genome\t" . $detected_genome/$total_genome . "\n";
	
	print STDERR "Overlap:  Genome\t$ov_len_genome\t$total_len_genome\t". $ov_len_genome/$total_len_genome . "\n\n";
}


exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "172.17.141.2";
        my $username = "sossowski";
        my $password = "";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

