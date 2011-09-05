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
#  Module: Analysis::Validation::PCR::ma_validation_ssr_variation.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $dbh;
&connect_to_db();

foreach my $file ("MA-30-119.SSR_variation.txt", "MA-30-29.SSR_variation.txt", "MA-30-49.SSR_variation.txt", "MA-30-59.SSR_variation.txt", "MA-30-69.SSR_variation.txt") {

open FILE, $file;


while( <FILE> ) {
	my @a = split " ", $_;
	my $chr    = $a[0];
	my $pos    = $a[1];
	my $pos2   = $a[2];
	my $start  = $pos - 500;
	my $end    = $pos2 + 500;
	my $seq    = $a[5];
	my @b = split '\.', $file;
	my $ma = $b[0];

	# Get reference sequence
	my $q2="SELECT	position, base 
		FROM 	seq_ref
		WHERE	chromosome = $chr && position BETWEEN $start AND $end
		ORDER BY chromosome, position\n";
	my $sth2 = $dbh->prepare($q2);
	$sth2->execute();

	my %cons_seq = ();
	my %good_ref = ();
	while( my $ref2 = $sth2->fetchrow_hashref() ) {
		$cons_seq{$ref2->{position}} = $ref2->{base};
		$good_ref{$ref2->{position}} = 0;
	}

	### Check reference call for all lines
	foreach my $table ("poly_reference_30_119", "poly_reference_30_29", "poly_reference_30_49",
				"poly_reference_30_59", "poly_reference_30_69") {

		### Get consensus sequence of surrounding region
		my $q3="SELECT 	position
			FROM 	$table
			WHERE	chromosome = $chr &&
				position BETWEEN $start AND $end
			ORDER by chromosome, position
		";
		my $sth3 = $dbh->prepare($q3);
		$sth3->execute();

		while( my $ref3 = $sth3->fetchrow_hashref() ) {
			$good_ref{$ref3->{position}}++;
		}
	}

	# Print sequence for primer design
	print ">$chr-$start-$end(SSR=MA_$ma-$pos-$pos2-$seq):TARGET=450,101\n";

	for(my $i = $start; $i <= $end; $i++) {
		if( (exists$good_ref{$i}) && ($good_ref{$i} == 5) ) {
			print $cons_seq{$i};
		}
		else { 
			if ($i >= $pos && $i <= $pos2) {
				print $cons_seq{$i};
			}
			else {
				print "N"; 
				#print $cons_seq{$i};
			}
		}
	}
	print "\n\n";
}

close FILE;

}

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ma";
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

