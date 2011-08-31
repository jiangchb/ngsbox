#!/usr/bin/perl

use strict;
use warnings;
use DBI;
use Getopt::Long;

my $database;
my $summary;
my $gc_low;
my $gc_up;
my $gc_mid_low;
my $gc_mid_up;
my $chr;

my %CMD;
GetCom();

my $dbh;
&connect_to_db();

my @COV_ALL = ();
my @COV_UP = ();
my @COV_LOW = ();
my @COV_MID = ();

my $q2 = "SELECT coverage, gc_cont FROM $summary where seg_flag = 'u' and ref_base <> 'N' and chromosome in ($chr)";
my $sth2 = $dbh->prepare($q2); $sth2->execute();

while (my $res = $sth2->fetchrow_hashref()) {
	if (defined($res->{gc_cont})) {
		if ($res->{gc_cont} <= $gc_low) {
			$COV_LOW[$res->{coverage}]++;
		}
		if ($res->{gc_cont} >= $gc_up) {
                	$COV_UP[$res->{coverage}]++;
	        }
		if ($res->{gc_cont} >= $gc_mid_low and $res->{gc_cont} <= $gc_mid_up) {
			$COV_MID[$res->{coverage}]++;
		}
		$COV_ALL[$res->{coverage}]++;
	}
}

open ALL, "> cov_hist.all.out" or die "Cannot open all.cov_hist.out\n";
open MID, "> cov_hist.mid.out" or die "Cannot open mid.cov_hist.out\n";
open LOW, "> cov_hist.low.out" or die "Cannot open low.cov_hist.out\n";
open UP, "> cov_hist.up.out" or die "Cannot open up.cov_hist.out\n";

for (my $i=0; $i<@COV_ALL; $i++) {
	if (defined($COV_ALL[$i])) {
		print ALL $i, "\t", $COV_ALL[$i], "\n";
	} else {
		print ALL $i, "\t0\n";
	}
}
for (my $i=0; $i<@COV_LOW; $i++) {
        if (defined($COV_LOW[$i])) {
                print LOW $i, "\t", $COV_LOW[$i], "\n";
        } else {
                print LOW $i, "\t0\n";
        }
}
for (my $i=0; $i<@COV_UP; $i++) {
        if (defined($COV_UP[$i])) {
                print UP $i, "\t", $COV_UP[$i], "\n";
        } else {
                print UP $i, "\t0\n";
        }
}
for (my $i=0; $i<@COV_MID; $i++) {
	if (defined($COV_MID[$i])) {
		print MID $i, "\t", $COV_MID[$i], "\n";
	}
	else {
		print MID $i, "\t0\n";
	}
}


close ALL or die "Cannot close all.cov_hist.out\n";
close LOW or die "Cannot close low.cov_hist.out\n";
close UP or die "Cannot close up.cov_hist.out\n";
close MID or die "Cannot close open up.cov_hist.out\n";


#system("R --slave --vanilla --args all.cov_hist.out low.cov_hist.out up.cov_hist.out < /ebio/abt6/korbinian/pgsp/Analysis/Panel_01/Fig01/coverage_hist.R");


exit(0);

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

sub GetCom {

  my @usage = ("\nUsage: $0 --consensus=<table> --ref=<table> --chr=<ids>

required:
--database\tDatabase name
--summary\tSummary table
--gc_low\tGC content defining the upper border for the low gc content set
--gc_up\tGC content defining the lower border for the high gc content set
--gc_mid_low\tGC content defining the lowerd border for the mid range gc content set
--gc_mid_up\tGC content defining the upper border for the mid range gc content set
--chr\tChromosomes, comma separated

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "summary=s","gc_low=s", "gc_up=s", "chr=s", "gc_mid_low=s", "gc_mid_up=s");

	die("Please specify database name\n") unless $CMD{database};
        die("Please specify summary table name\n") unless $CMD{summary};
        die("Please specify gc low threshold\n") unless $CMD{gc_low};
	die("Please specify gc up threshold\n") unless $CMD{gc_up};
	die("Please specify gc up threshold\n") unless $CMD{gc_mid_up};
	die("Please specify gc up threshold\n") unless $CMD{gc_mid_low};
	die("Please specify chromosomes\n") unless $CMD{chr};


	$database = $CMD{database};
	$summary = $CMD{summary};
	$gc_low = $CMD{gc_low};
	$gc_up = $CMD{gc_up};
	$gc_mid_low = $CMD{gc_mid_low};
	$gc_mid_up = $CMD{gc_mid_up};
	$chr = $CMD{chr};
}


