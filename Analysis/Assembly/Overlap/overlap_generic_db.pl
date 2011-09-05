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
#  Module: Analysis::Assembly::Overlap::overlap_generic_db.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

### User params
my $usage = "\n\n$0 ecotype number_chromosomes variants_min_length variants_max_length extension database table1 table2\n\n";
my $ecotype    = shift or die $usage;
my $num_chr    = shift or die $usage;
my $min_length = shift or die $usage;
my $max_length = shift or die $usage;
my $extension  = shift; #or die $usage;
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


for(my $chr = 1; $chr <= $num_chr; $chr++) {

	### Load ranged features of table 1
	my $q ="SELECT 	distinct chromosome, begin, end 
		FROM 	$table1
		WHERE 	ecotype = '$ecotype'
			&& chromosome = $chr
			# && quality >= 32
			# && (end - begin + 1) between $min_length AND $max_length
			# && len between $min_length AND $max_length
		ORDER BY chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %variants1 = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin} - $extension; $i <= $ref->{end} + $extension; $i++) {
			$variants1{$i} = 1;
		}
	}


	### Load ranged features of table 2 and compare
	$q = "	SELECT 	distinct chromosome, begin, end
		FROM	$table2
		WHERE 	ecotype = '$ecotype'
			&& chromosome = $chr
			&& len between $min_length AND $max_length
			# && (end - begin + 1) between $min_length AND $max_length
			# && support >= 12 && concordance >= 0.8 && repetitiveness = 1
			# && quality >= 25
		ORDER by chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		
		$total_chr{$chr}++;
		$total_genome++;

		$total_len_chr{$chr} += ($ref->{end} - $ref->{begin} + 1);
		$total_len_genome += ($ref->{end} - $ref->{begin} + 1);

		my $found = 0;
		for(my $pos = $ref->{begin}; $pos <= $ref->{end}; $pos++) {
			if(exists $variants1{$pos}) {
				$ov_len_chr{$chr}++;
				$ov_len_genome++;
				$found = 1;
			}
		}
		if( $found == 1) {
			$detected_chr{$chr}++;
			$detected_genome++;
		}
	}

	if( ($total_chr{$chr} > 0) && ($total_len_chr{$chr} > 0) ) {
		print "Detected: Chr$chr\t" . 
			$detected_chr{$chr} ."\t". $total_chr{$chr} ."\t". 
			$detected_chr{$chr}/$total_chr{$chr} . "\n";
		
		print "Overlap:  Chr$chr\t" . 
			$ov_len_chr{$chr} ."\t". $total_len_chr{$chr} ."\t". 
			$ov_len_chr{$chr}/$total_len_chr{$chr} . "\n\n";
	}
}

if( ($total_genome > 0) && ($total_len_genome > 0) ) {
	print "Detected: Genome\t$detected_genome\t$total_genome\t" . $detected_genome/$total_genome . "\n";
	
	print "Overlap:  Genome\t$ov_len_genome\t$total_len_genome\t". $ov_len_genome/$total_len_genome . "\n\n";
}


exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

