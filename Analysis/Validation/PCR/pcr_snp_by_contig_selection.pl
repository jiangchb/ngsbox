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
#  Module: Analysis::Validation::PCR::pcr_snp_by_contig_selection.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $database = shift;
my $file     = shift;
my $dbh;
&connect_to_db();
open IN, $file or die "Cannot open $file\n";

while(<IN>) {
	chomp;
	my ($eco, $contig, $chr, $begin, $end) = split(/\t/, $_);

	my $q ="SELECT chromosome, position 
		FROM poly_snp_hq 
		WHERE	ecotype = '$eco' &&
			contig = '$contig' &&
			chromosome = $chr &&
			position between $begin AND $end
		ORDER BY chromosome, position";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	print "$contig\t$chr\t$begin\t$end\tSNP: ";
	while( my $ref = $sth->fetchrow_hashref() ) {
		my $pos = $ref->{position};
		print "$pos,";
	}
	print "\n";
}

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

