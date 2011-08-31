#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use DBI;

my %CMD = ();
my $database;
my $samples;
my $file;

GetCom();

my @sample_array = split(/,/, $samples);
my $dbh;
&connect_to_db();

open IN, $file or die "Cannot open input file\n";
while( <IN> ) {
	chomp;

	my @entries = split(/\t/, $_);
	my $cov_fwd = 0;
	my $cov_rev = 0;

	foreach my $sample ( @sample_array ) {
		my $q ="SELECT 	coverage_fwd, coverage_rev
			FROM 	$sample
			WHERE 	length between 19 AND 22 &&
				chromosome = $entries[0] &&
				position between $entries[1] AND $entries[1] + 2
			";
		my $sth = $dbh->prepare($q); $sth->execute();
		my $ref = $sth->fetchrow_hashref();

		$cov_fwd += $ref->{coverage_fwd} if defined $ref->{coverage_fwd};
		$cov_rev += $ref->{coverage_rev} if defined $ref->{coverage_rev};
	}

	if($cov_fwd >= 1 || $cov_rev >= 1) {
		#print substr($entries[4], 0, 10) . "\n";
		print "$entries[4]\t$entries[0]\t$entries[1]\t$entries[2]\t$entries[3]\t".
			"$entries[6]\t$entries[5]\t$cov_fwd\t$cov_rev\n";
	}
}

exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\tDB name
--samples\t\tConsensus tables, comma delimeted
--file\t\tinput file with miRNA predictions

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "samples=s", "file=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify samples\n") unless defined($CMD{samples});
	die("Please specify file\n") unless defined($CMD{file});

	$database     = $CMD{database};
	$samples      = $CMD{samples};
	$file         = $CMD{file};

	return(0);
}

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


exit(0);
