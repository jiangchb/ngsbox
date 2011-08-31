#!/usr/bin/perl
use strict;
use warnings;
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

