#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;

my $database;
my $summary;

my %counts = ();

my %CMD;
GetCom();

my $dbh;
&connect_to_db();

open ALL, ">cov_hist.all.out" or die "Cannot open cov_hist.all.out\n";

my $q0 = "SELECT chromosome FROM seq_max WHERE chromosome BETWEEN 3 and 8";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();

while (my $ref0 = $sth0->fetchrow_hashref()) {
	my $chr = $ref0->{chromosome};

	my $q = "SELECT position, nonrep_coverage FROM $summary where chromosome = $chr";
	my $sth = $dbh->prepare($q); $sth->execute();

	my $last_pos = 0;
	while (my $ref = $sth->fetchrow_hashref()) {
		my $pos = $ref->{position};
		my $cov = $ref->{nonrep_coverage};
	
		if($pos != $last_pos + 1) {
			for(my $i = $last_pos + 1; $i < $pos; $i++) {
				$counts{0}++;
			}
		}
		$counts{$cov}++;
		$last_pos = $pos;
	}
}

foreach (sort keys %counts) {
	print ALL "$_\t$counts{$_}\n";
}

close ALL or die "Cannot close cov_hist.all.out\n";

system("R --slave --vanilla < /ebio/abt6/stephan/pgsp/ChipSeq/coverage_hist.R");

exit(0);



### Connect to mysql database ###
sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}



### Get comman line parameters ###
sub GetCom {

  my @usage = ("\nUsage: $0 --database=<table> --summary=<table> --chr=<ids>

required:
--database\tDatabase name
--summary\tSummary table

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "summary=s");

	die("Please specify database name\n") unless $CMD{database};
        die("Please specify summary table name\n") unless $CMD{summary};

	$database = $CMD{database};
	$summary = $CMD{summary};
}


