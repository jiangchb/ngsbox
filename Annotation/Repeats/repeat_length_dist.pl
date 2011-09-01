#!/usr/bin/perl -w

use strict;
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
