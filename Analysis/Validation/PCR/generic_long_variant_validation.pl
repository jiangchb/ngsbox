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
#  Module: Analysis::Validation::PCR::generic_long_variant_validation.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $file = shift;
my $callability_check = shift;
my $target_Xed = shift;
my $extend = shift;

my $dbh;
&connect_to_db();

open FILE, $file or die "\n\nCannot open infile\n\n";
while( <FILE> ) {
	my @a = split " ", $_;
	my $strain         = $a[0];
	my $chr            = $a[1];
	my $beg_variant    = $a[2];
	my $end_variant    = $a[3];
	my $variant        = $a[4];

	my $beg_sanger  = $beg_variant - $extend;
	my $end_sanger  = $end_variant + $extend;


	### Get reference sequence
	my $q= "SELECT	position, base
		FROM 	seq_ref
		WHERE	chromosome = $chr && position BETWEEN $beg_sanger AND $end_sanger
		ORDER BY chromosome, position\n";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %cons_seq = ();
	my %good_ref = ();
	while( my $ref = $sth->fetchrow_hashref() ) {
		$cons_seq{$ref->{position}} = $ref->{base};
		$good_ref{$ref->{position}} = 0;
	}


	### Check reference call for all lines
	foreach my $table ("poly_reference_30_119", "poly_reference_30_29", "poly_reference_30_49", "poly_reference_30_59", "poly_reference_30_69") {

		### Get consensus sequence of surrounding region
		my $q2="SELECT 	position
			FROM 	$table
			WHERE	chromosome = $chr &&
				position BETWEEN $beg_sanger AND $end_sanger
			ORDER by chromosome, position
		";
		my $sth2 = $dbh->prepare($q2);
		$sth2->execute();

		while( my $ref2 = $sth2->fetchrow_hashref() ) {
			$good_ref{$ref2->{position}}++;
		}
	}


	### Print header
	my $target_start = $beg_variant - $beg_sanger;
	my $target_length = $end_variant - $beg_variant + 1;
	print ">$chr-$beg_sanger-$end_sanger(VAR=$strain-$beg_variant-$end_variant-$variant):TARGET=$target_start,$target_length\n";

	### Print sequence
	for(my $i = $beg_sanger; $i <= $end_sanger; $i++) {
		if( $i >= $beg_variant && $i <= $end_variant && $target_Xed == 1 ) {
			print "X";
		}
		elsif( ($callability_check == 0) || (exists$good_ref{$i} && $good_ref{$i} == 5) ) {
			print $cons_seq{$i};
		}
		else { 
			print "N"; 
		}
	}
	print "\n\n";
}

close FILE;

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ma";
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

