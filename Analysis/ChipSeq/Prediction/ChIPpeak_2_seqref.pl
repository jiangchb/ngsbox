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
#  Module: Analysis::ChipSeq::Prediction::ChIPpeak_2_seqref.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $type = shift;
my $dbh;
&connect_to_db();

### Query chromosome length
my $q = "SELECT chromosome, chr_alias FROM seq_max where chromosome between 3 and 8 ORDER BY chromosome";
my $sth = $dbh->prepare($q);
$sth->execute();

### Parse TFBS predictions 
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr       = $ref->{chromosome};
        my $chr_alias = $ref->{chr_alias};
	my %tfbs = ();

	### Load TFBS prediction
	foreach my $tfbs_table ("binding_site", "cluster", "redfly_tfbs", "redfly_crm") {
		my $q= "SELECT 	begin, end
			FROM 	$tfbs_table
			WHERE 	chromosome = '$chr_alias'
			ORDER BY begin
		";
		my $sth = $dbh->prepare($q);
		$sth->execute();
	
		while(my $ref = $sth->fetchrow_hashref()) {
			for(my $i = $ref->{begin} - 100; $i <= $ref->{end} + 100; $i++) {
				$tfbs{$i} = 1;
			}
		}
	}

	### Load ChipSeq peaks
	my $q="	SELECT 	distinct begin, end
		FROM	enriched_dfd_s1_nonrep_gene
		WHERE 	chromosome = $chr &&
			max_cov >= 6 &&
			segment_length >= 75 &&
			fbgn in (select sequence_derived_from from array_dfd_anno_gene_jun08 where dfd_early_no_rec_hs >= 1.5 || dfd_early_no_rec_hs <= 0.4 || dfd_early_gof > 1.5 || dfd_early_gof between 0.001 and 0.5)
		ORDER by chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};
		my $hit = 0;

		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $tfbs{$pos}) { $hit = 1; }
		}
		#if( ($hit == 1) && ($type == 1) ) { 
			&print_seqref($chr, $chr_alias, $beg, $end, "dfd");
		#}
		#elsif( ($hit == 0) && ($type == 0) ) {
			&print_seqref($chr, $chr_alias, $beg, $end, "nodfd");
		#}
	}
}

exit(0);

sub print_seqref
{
	my ($chr, $chr_alias, $beg, $end, $name) = @_;
	my $q="SELECT base FROM seq_ref WHERE chromosome = $chr && position BETWEEN $beg - 100 AND $end - 100 ORDER BY position";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	print ">$name"."_loc=$chr_alias:$beg..$end\n";
	while(my $ref = $sth->fetchrow_hashref()) {
		print $ref->{base};
	}
	print "\n";
}

#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "chip_seq_dmel";
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

