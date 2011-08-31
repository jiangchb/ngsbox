#!/usr/bin/perl

use strict;
use warnings;
use DBI;

my $ecotype = shift;

### Query fragments from 2010 set
my $dbh;
&connect_to_db();
my $q = "SELECT chromosome, position,end FROM mn_fragments
		WHERE ecotype = '$ecotype'
		ORDER BY chromosome, position";
my $sth = $dbh->prepare($q);
$sth->execute();

while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr = $ref->{chromosome};
	my $pos = $ref->{position};
	my $end = $ref->{end};

	my $q2 = "SELECT position, seg_flag FROM consensus_bur 
		WHERE chromosome = $chr &&
			position BETWEEN $pos AND $end
		ORDER BY position";
	my $sth2 = $dbh->prepare($q2);
	$sth2->execute();

	while( my $ref2 = $sth2->fetchrow_hashref() ) {
		print "$ecotype\t$chr\t$ref2->{position}\t$ref2->{seg_flag}\n";
	}
}


$dbh->disconnect();
exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
	my $databaseName = "solexa";
	my $driver = "mysql";
	my $host = "ume.fml.local";
	my $username = "solexa";
	my $password = "s0lexa";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}
