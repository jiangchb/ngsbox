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
#  Module: Analysis::SV::SV_Validation::overlap_Orphans_Unseq.pl
#  Purpose:
#  In:
#  Out:
#


use DBI;

my $ecotype    = shift;

### Connect to database
my $dbh;
&connect_to_db();


for(my $chr = 1; $chr <= 1; $chr++) {

	### Load Unsequenced regions
	my $q ="SELECT 	begin, end
		FROM 	unsequenced
		WHERE 	sample = '$ecotype' &&
			chromosome = $chr
		ORDER BY chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %unseq = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		$unseq{$ref->{begin}} = $ref->{end};
	}


        ### Load core unsequenced regions
        $q = "  SELECT  begin, end
                FROM    unsequenced_core
                WHERE   sample = '$ecotype' &&
                        chromosome = $chr
                ORDER BY chromosome, begin
        ";
        $sth = $dbh->prepare($q);
        $sth->execute();

        my %unseq_core = ();
        while(my $ref = $sth->fetchrow_hashref()) {
		$unseq_core{$ref->{begin}} = $ref->{end};
        }

        ### Load core nonrep unsequenced regions
        $q = "  SELECT  begin, end
                FROM    unsequenced_cn
                WHERE   sample = '$ecotype' &&
                        chromosome = $chr
                ORDER BY chromosome, begin
        ";
        $sth = $dbh->prepare($q);
        $sth->execute();

        my %unseq_cn = ();
        while(my $ref = $sth->fetchrow_hashref()) {
		$unseq_cn{$ref->{begin}} = $ref->{end};
        }


	### Load orphan reads
	$q = "  SELECT  begin, end, pattern_support
		FROM    sv_missing_ends
		WHERE   sample = '$ecotype' &&
			chromosome = $chr
		ORDER by chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %orphans = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$orphans{$i} = $ref->{pattern_support};
		}
	}


	### Compare
	foreach my $unseq_beg (sort {$a<=>$b} keys %unseq_cn) {
		my $unseq_end = $unseq_cn{$unseq_beg};
		my $unseq_len = $unseq_end - $unseq_beg + 1;

		my $orphan_support = 0;
		for(my $pos = $unseq_beg; $pos <= $unseq_end; $pos++) {
			if(exists $orphans{$pos}) {
				$orphan_support = $orphans{$pos};
				last;
			}
		}

		if($orphan_support > 0.1) {
			print "$ecotype\t$chr\t$unseq_beg\t$unseq_end\t$unseq_len\t$orphan_support\n";
		}
	}
}
exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ath_pe";
        my $driver = "mysql";
        my $host = "ume.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

