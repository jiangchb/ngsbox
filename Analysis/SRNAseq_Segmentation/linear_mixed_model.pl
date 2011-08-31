#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use DBI;

my %CMD = ();
my $database;
my $consensus_tables;
my $sample_names;
my $srna_species;
my $cov_type;
my $win_size;

GetCom();
my $dbh;
&connect_to_db();

my @consensus_tables = split(",", $consensus_tables);
my @sample_names = split(/,/, $sample_names);
my %sample_desc = ();

### Parse sample description
my $q = "SELECT table_name,sample,flowcell,lane,reads_24,assay,organ,biorep,techrep,temperature FROM library_translation";
my $sth = $dbh->prepare($q);
$sth->execute();
while( my $ref = $sth->fetchrow_hashref() ) {
	my $table_name = $ref->{table_name};
	my $entries = $ref->{sample} ."\t". $ref->{flowcell} ."\t". $ref->{lane} ."\t". $ref->{reads_24} ."\t". $ref->{assay} ."\t". $ref->{organ} ."\t". $ref->{biorep} ."\t". $ref->{techrep} ."\t". $ref->{temperature};
	$sample_desc{$table_name} = $entries;
}


### Query chromosome length
$q = "SELECT chromosome, max_pos FROM seq_max ORDER BY chromosome";
$sth = $dbh->prepare($q);
$sth->execute();

### Parse each chromosome
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr = $ref->{chromosome};
	my $chr_size = $ref->{max_pos};
	my  $sum_cov_sample1 = 0;
	my  $sum_cov_sample2 = 0;

	### Parse each bin in chromosome
	for (my $window = 1; $window <= $chr_size; $window+=$win_size) {
		my @sample_cov = ();

		### Parse each sample
		for ( my $i = 0; $i < @consensus_tables; $i++ ) {
			$sample_cov[$i] = 0;

			my $q2="SELECT 	$cov_type, average_hits
				FROM 	$consensus_tables[$i]
        	        	WHERE 	chromosome = $chr &&
					position BETWEEN $window AND ($window + $win_size - 1) &&
					length = $srna_species
	        		";
			my $sth2 = $dbh->prepare($q2); $sth2->execute();
			while (my $res = $sth2->fetchrow_hashref()) {
				$sample_cov[$i]  += ( $res->{$cov_type} / $res->{average_hits} );
			}

			### Print correlation plot
			print $sample_desc{$consensus_tables[$i]} . "\t$sample_cov[$i]\n";
		}
	}
}

exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\t\tDB name
--consensus_tables\t\tConsensus table names, comma seperated (no spaces allowed)
--sample_names\t\tName of all samples, comma seperated (no spaces allowed)
--srna_species\t\tsRNA length
--cov_type\t\tCoverage type used to segment (coverage, coverage_forward ...)
--win_size\t\tSize of sliding window

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "consensus_tables=s", "sample_names=s", "srna_species=s", "cov_type=s", "win_size=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify consensus_tables\n") unless defined($CMD{consensus_tables});
	die("Please specify sample_names\n") unless defined($CMD{sample_names});
	die("Please specify srna_species\n") unless defined($CMD{srna_species});
	die("Please specify cov_type\n") unless defined($CMD{cov_type});
	die("Please specify win_size\n") unless defined($CMD{win_size});

	$database         = $CMD{database};
	$consensus_tables = $CMD{consensus_tables};
	$sample_names     = $CMD{sample_names};
	$srna_species     = $CMD{srna_species};
	$cov_type         = $CMD{cov_type};
	$win_size         = $CMD{win_size};

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
