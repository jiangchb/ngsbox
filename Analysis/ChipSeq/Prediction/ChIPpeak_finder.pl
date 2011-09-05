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
#  Module: Analysis::ChipSeq::Prediction::ChIPpeak_finder.pl
#  Purpose:
#  In:
#  Out:
#

use Getopt::Long;
use DBI;

my %CMD = ();
my $database;
my $sample;
my $win_size;
my $cov_type;
my $cov_min;
my $peak_min;
my $min_len;

GetCom();

my $dbh;
&connect_to_db();

### Query chromosome length
my $q = "SELECT chromosome, chr_alias FROM seq_max ORDER BY chromosome";
my $sth = $dbh->prepare($q);
$sth->execute();

### Parse Segments
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr       = $ref->{chromosome};
	my $chr_alias = $ref->{chr_alias};
	
	my $start     = 0;
	my $last_pos  = 0;
	my $sum_cov   = 0;
	my $max_cov   = 0;

	my $q2="SELECT 	position, $cov_type
		FROM 	$sample
                WHERE 	chromosome = $chr AND 
			$cov_type >= $cov_min
	";
	my $sth2 = $dbh->prepare($q2); $sth2->execute();
	while (my $res = $sth2->fetchrow_hashref()) {
		my $pos = $res->{position};

		### Print segment
		if( $pos > ($last_pos + $win_size) ) {
			my $seg_length = $last_pos - $start + 1;
			my $avg_cov = $sum_cov / $seg_length;

			if( ($max_cov >= $peak_min) && ($seg_length > $min_len) ) {
				print "$chr\t$chr_alias\t$start\t$last_pos\t$seg_length\t" .
					sprintf("%.2f", $avg_cov) .
					"\t$max_cov\n";
			}
			$start   = $pos;
			$sum_cov = 0;
			$max_cov = 0;
		}

		### Update segment features
		$last_pos = $pos;
		$sum_cov += $res->{$cov_type};
		if( $res->{$cov_type} > $max_cov ) {
			$max_cov = $res->{$cov_type};
		}
	}
}

exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\tDB name
--sample\t\tConsensus table name
--win_size\t\tSize of sliding window
--cov_type\t\tcoverage, coverage_forward, coverage_reverse, nonrep_...
--cov_min\t\tCoverage threshold
--peak_min\t\tMinimum peak coverage
--min_len\t\tMinimum segment length
\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "sample=s", "win_size=s", "cov_type=s", "cov_min=s", "peak_min=s", "min_len=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify sample\n") unless defined($CMD{sample});
	die("Please specify win_size\n") unless defined($CMD{win_size});
	die("Please specify cov_type\n") unless defined($CMD{cov_type});
	die("Please specify cov_min\n") unless defined($CMD{cov_min});
	die("Please specify peak_min\n") unless defined($CMD{peak_min});
	die("Please specify min_len\n") unless defined($CMD{min_len});

	$database        = $CMD{database};
	$sample          = $CMD{sample};
	$win_size        = $CMD{win_size};
	$cov_type        = $CMD{cov_type};
	$cov_min         = $CMD{cov_min};
	$peak_min        = $CMD{peak_min};
	$min_len         = $CMD{min_len};

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
