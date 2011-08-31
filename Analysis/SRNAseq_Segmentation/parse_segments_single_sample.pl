#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use DBI;

my %CMD = ();
my $database;
my $sample;
my $srna_length_min;
my $srna_length_max;
my $cov_min;
my $win_size;

GetCom();

my $dbh;
&connect_to_db();

foreach my $cov_type ("coverage_fwd", "coverage_rev") {
	for (my $length = $srna_length_min; $length <= $srna_length_max; $length++) {
		### Query chromosome length
		my $q = "SELECT chromosome, max_pos FROM seq_max ORDER BY chromosome";
		my $sth = $dbh->prepare($q);
		$sth->execute();

		### Parse Segments
		while( my $ref = $sth->fetchrow_hashref() ) {
			my $chr = $ref->{chromosome};
			my $chr_size = $ref->{max_pos};
	
			my $start       = 0;
			my $last_pos    = 0;
			my $sum_read    = 0;
			my $sum_avg_hit = 0;
			my $max_read_start = 0;

			my $q2="SELECT 	position, $cov_type, average_hits
				FROM 	$sample
		                WHERE 	chromosome = $chr AND 
					length = $length &&
					$cov_type > 0
			";
			my $sth2 = $dbh->prepare($q2); $sth2->execute();
			while ( my $res = $sth2->fetchrow_hashref() ) {
				my $pos = $res->{position};
	
				if( $pos > ($last_pos + $win_size) ) {
					if($sum_read >= $cov_min) {
						$last_pos += $length;
						my $seg_length = $last_pos - $start;
						my $avg_cov = $sum_read / $seg_length * $length;
						my $avg_avg_hit = $sum_avg_hit / $sum_read;
						my $orientation = "+";
						if($cov_type eq "coverage_rev") { $orientation = "-"; }

						print 	"$length\t$orientation\t$chr\t$start\t$last_pos\t".
							"$seg_length\t$sum_read\t$max_read_start\t".
							sprintf("%.2f", $avg_cov) ."\t". 
							sprintf("%.2f", $avg_avg_hit) ."\n";
					}
					$start = $pos;
					$sum_read = 0;
					$sum_avg_hit = 0;
					$max_read_start = 0;
				}
				$last_pos = $pos;
				$sum_read += $res->{$cov_type};
				$sum_avg_hit += ($res->{average_hits} * $res->{$cov_type});
				if($res->{$cov_type} > $max_read_start) { $max_read_start = $res->{$cov_type}; }
			}
		}
	}
}

exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\tDB name
--sample\t\tConsensus table name
--srna_length_min\t\tMin sRNA species length
--srna_length_max\t\tMax sRNA species length
--cov_min\t\tMin. coverage to report position
--win_size\t\tSize of sliding window

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "sample=s", "srna_length_min=s", "srna_length_max=s", "cov_min=s", "win_size=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify sample\n") unless defined($CMD{sample});
	die("Please specify srna_length_min\n") unless defined($CMD{srna_length_min});
	die("Please specify srna_length_max\n") unless defined($CMD{srna_length_max});
	die("Please specify cov_min\n") unless defined($CMD{cov_min});
	die("Please specify win_size\n") unless defined($CMD{win_size});

	$database        = $CMD{database};
	$sample          = $CMD{sample};
	$srna_length_min = $CMD{srna_length_min};
	$srna_length_max = $CMD{srna_length_max};
	$cov_min         = $CMD{cov_min};
	$win_size        = $CMD{win_size};

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
