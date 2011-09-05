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
#  Module: Analysis::Validation::PCR::pcr_validation.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $database = shift;
my $ecotype  = shift;
my $type     = shift;
my $support  = shift;

my $dbh;
&connect_to_db();

if( $type eq "SNP" ) {
	### Get SNPs
	my $q ="SELECT contig, chromosome, position 
		FROM poly_snp_hq
		WHERE   ecotype = '$ecotype' &&
			support >= $support &&
			based_on = 'core nonrep'
		ORDER by chromosome, position";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	while( my $ref = $sth->fetchrow_hashref() ) {
		my $random_number = int(rand(5000));
		if($random_number < 1500) {
			my $contig = $ref->{contig};
			my $chr    = $ref->{chromosome};
			my $pos    = $ref->{position};
			my $start  = $pos - 500;
			my $end    = $pos + 500;
	
			### Get consensus sequence of surrounding region
			my $q2="SELECT 	position, base_call, repeat_type
				FROM 	poly_reference
				WHERE	chromosome = $chr &&
					position BETWEEN $start AND $end
				ORDER by chromosome, position
			";
			my $sth2 = $dbh->prepare($q2);
			$sth2->execute();
	
			my %cons_seq = ();
			my $no_ambi_counter = 0;
			my $no_repeat_counter = 0;
			while( my $ref2 = $sth2->fetchrow_hashref() ) {
				$cons_seq{$ref2->{position}} = $ref2->{base_call};
				$no_ambi_counter++;
				if($ref2->{repeat_type} eq "u") { $no_repeat_counter++; }
			}
			if ($no_ambi_counter > 950 && $no_repeat_counter > 900) {
				print ">$contig-$chr-$start-$end($pos):TARGET=550,100\n";
				for( my $i = $start; $i <= $end; $i++) {
					if( exists $cons_seq{$i} ) { print $cons_seq{$i}; }
					else { print "N"; }
				}
				print "\n";
			}
		}
	}
}

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

