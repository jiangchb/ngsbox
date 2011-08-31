#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $dbh;
&connect_to_db();

my $q0 = "SELECT chromosome FROM seq_max ORDER BY chromosome";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();

while(my $ref0 = $sth0->fetchrow_hashref()) {
	my $chr = $ref0->{chromosome};
	print ">$chr\n";

	my $q = "SELECT GC_content_101 FROM seq_ref WHERE chromosome = $chr ORDER BY position";
	my $sth = $dbh->prepare($q);
	$sth->execute();
	while(my $ref = $sth->fetchrow_hashref()) {
		if($ref->{GC_content_101}) {
			print chr(int($ref->{GC_content_101}) + 34);
		}
		else {
			print chr(33);
		}
	}
	$sth->finish();

	print "\n";
}
$sth0->finish();

exit(0);

##############################################################################
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
		
