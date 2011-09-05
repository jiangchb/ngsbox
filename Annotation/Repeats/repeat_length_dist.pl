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
#  Module: Annotation::Repeats::repeat_length_dist.pl
#  Purpose:
#  In:
#  Out:
#


use DBI;

my $usage = "repeat_length_dist.pl\n";
my $database = "solexa";

my $dbh;
&connect_to_db();

my $q = "SELECT chromosome, position, seg_flag from consensus_col ORDER BY chromosome, position";
my $sth = $dbh->prepare($q);
$sth->execute();

my $curr_chr = -1;
my $rep_length = 0;
my $uniq_length = 0;

open REP, "> replength.txt";
open UNIQ, "> uniqlength.txt";

while(my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chromosome};
	my $pos = $ref->{position};
	my $seg = $ref->{seg_flag};

	if ($chr != $curr_chr) {
		if ($rep_length > 0) {
			print REP $rep_length, "\n";
		}
		if ($uniq_length > 0) {
			print UNIQ $uniq_length, "\n";
		}
	}
	$curr_chr = $chr;

	if ($seg eq "u") {
		if ($rep_length > 0) {
			print REP $rep_length, "\n";
		}
		$rep_length = 0;
		$uniq_length++;
	}
	else {
		if ($uniq_length > 0) {
                        print UNIQ $uniq_length, "\n";
                }
                $uniq_length = 0;
		$rep_length++;
	}
}
if ($rep_length > 0) {
	print REP $rep_length, "\n";
}
if ($uniq_length > 0) {
	print UNIQ $uniq_length, "\n";
}




#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db {
	my $databaseName = $database;
	my $driver = "mysql";
	my $host = "ume.eb.local";
	my $username = "solexa";
	my $password = "s0lexa";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh = DBI->connect($dsn, $username, $password ) || db_connect_error();
}
