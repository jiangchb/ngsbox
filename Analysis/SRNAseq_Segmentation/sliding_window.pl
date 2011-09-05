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
#  Module: Analysis::SRNAseq_Segmentation::sliding_window.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;
use DBI;

my $dbh;

my %CMD = ();
my $database;
my $sample1;
my $sample2;
my $length_min;
my $length_max;
my $cov_type;
my $win_size;
my $sample1_name;
my $sample2_name;
my $prefix;


GetCom();
&connect_to_db();

for (my $length = $length_min; $length <= $length_max; $length++) {
	my $file_name = "$prefix.$sample1_name.$sample2_name.l$length.bin$win_size";
	open OUT, ">$file_name.txt" or die "Cannot open output file\n";

        ### Query chromosome length
        my $q = "SELECT chromosome, max_pos FROM seq_max ORDER BY chromosome";
        my $sth = $dbh->prepare($q);
        $sth->execute();

        ### Parse Segments
        while( my $ref = $sth->fetchrow_hashref() ) {
		my $chr = $ref->{chromosome};
		my $chr_size = $ref->{max_pos};
		my %sample1_hash = ();
		my %sample2_hash = ();
		
		my $q2="SELECT 	position, $cov_type
			FROM 	$sample1
        	       	WHERE 	chromosome = $chr &&
				length= $length
		";
		my $sth2 = $dbh->prepare($q2); $sth2->execute();
		while (my $res = $sth2->fetchrow_hashref()) {
			$sample1_hash{$res->{position}}  = $res->{$cov_type};
		}


		$q2 = "	SELECT 	position, $cov_type
			FROM    $sample2
			WHERE 	chromosome = $chr &&
				length= $length
		";
                $sth2 = $dbh->prepare($q2); $sth2->execute();
                while (my $res = $sth2->fetchrow_hashref()) {
			$sample2_hash{$res->{position}}  += $res->{$cov_type};
                }


		my $sum_sample1 = 0;
		my $sum_sample2 = 0;
		for (my $win_pos = $win_size; $win_pos <= $chr_size - $win_size; $win_pos++) {
			my $sum_sample1 = 0;
			my $sum_sample2 = 0;

			for (my $pos = $win_pos - $win_size; $pos <= $win_pos + $win_size; $pos++) {
				$sum_sample1 += $sample1_hash{$pos} if defined $sample1_hash{$pos};
				$sum_sample2 += $sample2_hash{$pos} if defined $sample2_hash{$pos};
			}

			### Print correlation plot
			if($sum_sample1 == 0) { $sum_sample1 = 1;}
			if($sum_sample2 == 0) { $sum_sample2 = 1;}

			if($sum_sample1 > 1 || $sum_sample2 > 1) {
				print OUT "$sum_sample1\t$sum_sample2\n";
			}
		}
	}

	close OUT;
	system("sort -n -k1 -k2 $file_name.txt | uniq > $file_name.tmp");
	system("mv $file_name.tmp $file_name.txt");
	system("R --slave --vanilla --args $file_name $sample1_name $sample2_name < /ebio/abt6/stephan/pgsp/SRNA/Segmentation/coverage_correlation.R");
}


exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\t\tDB name
--sample1\t\tConsensus table name
--sample2\t\tConsensus table name
--length_min\t\tsRNA length min
--length_max\t\tsRNA length max
--cov_type\t\tCoverage type used to segment (coverage, coverage_forward ...)
--win_size\t\tSize of sliding window
--sample1_name\t\tName of sample 1 for output file name and plot description
--sample2_name\t\tName of sample 2 for output file name and plot description
--prefix\t\tPrefix for the output file name

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "sample1=s", "sample2=s", "length_min=s", "length_max=s", "cov_type=s", "win_size=s", "sample1_name=s", "sample2_name=s", "prefix=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify sample1\n") unless defined($CMD{sample1});
	die("Please specify sample2\n") unless defined($CMD{sample2});
	die("Please specify length_min\n") unless defined($CMD{length_min});
	die("Please specify length_max\n") unless defined($CMD{length_max});
	die("Please specify cov_type\n") unless defined($CMD{cov_type});
	die("Please specify win_size\n") unless defined($CMD{win_size});
	die("Please specify name of sample 1\n") unless defined($CMD{sample1_name});
	die("Please specify name of sample 2\n") unless defined($CMD{sample2_name});
	die("Please specify prefix\n") unless defined($CMD{prefix});

	$database     = $CMD{database};
	$sample1      = $CMD{sample1};
	$sample2      = $CMD{sample2};
	$length_min   = $CMD{length_min};
	$length_max   = $CMD{length_max};
	$cov_type     = $CMD{cov_type};
	$win_size     = $CMD{win_size};
	$sample1_name = $CMD{sample1_name};
	$sample2_name = $CMD{sample2_name};
	$prefix       = $CMD{prefix};

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
