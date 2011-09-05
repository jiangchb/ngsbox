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
#  Module: Analysis::SV::SV_Validation::overlap_SVins_Orphans_Unseq.pl
#  Purpose:
#  In:
#  Out:
#


use DBI;

my $ecotype    = shift;

### Connect to database
my $dbh;
&connect_to_db();


for(my $chr = 1; $chr <= 5; $chr++) {

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
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq{$i} = 1;
		}
	}

	### Load core unsequenced regions
	$q = "	SELECT  begin, end
		FROM    unsequenced_core
		WHERE   sample = '$ecotype' &&
			chromosome = $chr
		ORDER BY chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %unseq_core = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq_core{$i} = 1;
		}
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
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq_cn{$i} = 1;
		}
	}


	### Load SV insertion
	$q = "	SELECT 	begin, end, support, length, pvalue, spanned
		FROM	sv_insertion
		WHERE 	sample = '$ecotype' &&
			chromosome = $chr &&
			spanned > 1 && support < 40
		ORDER by chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %sv = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};
		my @features = ($ref->{end}, $ref->{support}, $ref->{'length'}, $ref->{pvalue}, $ref->{spanned});
		$sv{$beg} = \@features;
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
	my $total_support_cn = 0;
	my $total_count   = 0;
		
	foreach my $sv_beg (sort {$a<=>$b} keys %sv) {
		my $sv_end = $sv{$sv_beg}[0];
		$total_count++;

		# Unseq support for Insertion
		my $support_core = 0;
		my $support_cn = 0;
		for(my $pos = $sv_beg - 200; $pos <= $sv_end + 200; $pos++) {
			if(exists $unseq_core{$pos}) {
				$support_core = 1;
			}
			if(exists $unseq_cn{$pos}) {
				$support_cn = 1;
			}
		}

		# Orphan support for insertions
		my $orphan_support = 0;
		for(my $pos = $sv_beg; $pos <= $sv_end; $pos++) {
			if(exists $orphans{$pos}) {
				$orphan_support = $orphans{$pos};
				last;
			}
		}

		if(($support_cn == 1) || ($orphan_support > 0) ) { $total_support_cn++; }

		print "$ecotype\t$chr\t$sv_beg\t" . $sv{$sv_beg}[0] . "\t" . $sv{$sv_beg}[1] . "\t" . $sv{$sv_beg}[2] . "\t" . 
			$sv{$sv_beg}[3] . "\t" . $sv{$sv_beg}[4] . "\t$support_core\t$support_cn\t$orphan_support\n";
	}

	if($total_count != 0) {
		print STDERR "$chr\t$total_support_cn\t$total_count\t" . $total_support_cn/$total_count . "\n";
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

