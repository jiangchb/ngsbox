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
#  Module: Analysis::SRNAseq_Segmentation::parse_segments_multi_sample.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;
use DBI;

my $dbh;

my %CMD = ();
my $database;
my $coverage;;
my $sample1;
my $sample2;
my $length;

GetCom();
&connect_to_db();

my %COV_HASHES = ();

my %CHRSIZES = ();
$CHRSIZES{1} = 30432563;
$CHRSIZES{2} = 19705359;
$CHRSIZES{3} = 23470805;
$CHRSIZES{4} = 18585042;
$CHRSIZES{5} = 26992728;

for (my $chr = 1; $chr <= 5; $chr++) {

	my $chr_size = $CHRSIZES{$chr};
	
	for (my $batch = 1; $batch <= $chr_size; $batch+=100000) {

		my %SAMP1_FOR = ();
                my %SAMP1_REV = ();
		my %SAMP1_AVG_HITS = ();
                my %SAMP1_FOR_NONREP = ();
                my %SAMP1_REV_NONREP = ();
		my %SAMP1_AVG_HITS_NONREP = ();

		my %SAMP2_FOR = ();
                my %SAMP2_REV = ();
		my %SAMP2_AVG_HITS = ();
                my %SAMP2_FOR_NONREP = ();
                my %SAMP2_REV_NONREP = ();
		my %SAMP2_AVG_HITS_NONREP = ();

		my $q = "SELECT position, coverage_forward, coverage_reverse, average_hits, nonrep_coverage_forward, nonrep_coverage_reverse, nonrep_average_mismatches FROM $sample1
        	        	WHERE chromosome = $chr and position between $batch and $batch+100000-1 and coverage >= $coverage and length= $length
	        	        ";
		my $sth = $dbh->prepare($q); $sth->execute();
		while (my $res = $sth->fetchrow_hashref()) {
			$SAMP1_FOR{$res->{"position"}} = $res->{"coverage_forward"};
			$SAMP1_REV{$res->{"position"}} = $res->{"coverage_reverse"};
			$SAMP1_AVG_HITS{$res->{"position"}} = $res->{"average_hits"};
			$SAMP1_FOR_NONREP{$res->{"position"}} = $res->{"nonrep_coverage_forward"};
                        $SAMP1_REV_NONREP{$res->{"position"}} = $res->{"nonrep_coverage_reverse"};
			$SAMP1_AVG_HITS_NONREP{$res->{"position"}} = $res->{"nonrep_average_hits"};
		}


		$q = "SELECT position, coverage_forward, coverage_reverse, average_hits, nonrep_coverage_forward, nonrep_coverage_reverse, nonrep_average_mismatches FROM $sample2
                                WHERE chromosome = $chr and position between $batch and $batch+100000-1 and coverage >= $coverage and length = $length
                                ";
                $sth = $dbh->prepare($q); $sth->execute();
                while (my $res = $sth->fetchrow_hashref()) {
			$SAMP2_FOR{$res->{"position"}} = $res->{"coverage_forward"};
                        $SAMP2_REV{$res->{"position"}} = $res->{"coverage_reverse"};
			$SAMP2_AVG_HITS{$res->{"position"}} = $res->{"average_hits"};
                        $SAMP2_FOR_NONREP{$res->{"position"}} = $res->{"nonrep_coverage_forward"};
                        $SAMP2_REV_NONREP{$res->{"position"}} = $res->{"nonrep_coverage_reverse"};
			$SAMP2_AVG_HITS_NONREP{$res->{"position"}} = $res->{"nonrep_average_hits"};
                }
		
		for (my $i = $batch; $i<$batch+100000; $i++) {
			if (defined($SAMP1_FOR{$i}) or defined($SAMP2_FOR{$i})) {
				print $chr, "\t", $i, "\t";
				if (defined($SAMP1_FOR{$i})) { 
					print $SAMP1_FOR{$i}, "\t", $SAMP1_REV{$i}, "\t", $SAMP1_FOR_NONREP{$i}, "\t", $SAMP1_REV_NONREP{$i}, "\t", $SAMP1_AVG_HITS{$i}, "\t";
				}
				else {
					print "0\t0\t0\t0\t0\t";
				}
				if (defined($SAMP2_FOR{$i})) {
					print $SAMP2_FOR{$i}, "\t", $SAMP2_REV{$i}, "\t", $SAMP2_FOR_NONREP{$i}, "\t", $SAMP2_REV_NONREP{$i}, "\t", $SAMP2_AVG_HITS{$i}, "\n";
				}
				else {
					print "0\t0\t0\t0\t0\n";
				}
			}
		}
	}
}


exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\tDB name
--sample1\t\tConsensus table name
--sample2\t\tConsensus table name
--coverage\t\tMin. coverage in at least one sample to report position (5)
--length\t\tsRNA length

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "sample1=s", "sample2=s", "coverage=s", "length=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify sample1\n") unless defined($CMD{sample1});
	die("Please specify sample2\n") unless defined($CMD{sample2});
	die("Please specify coverage\n") unless defined($CMD{coverage});
	die("Please specify length\n") unless defined($CMD{length});

	$database = $CMD{database};
	$coverage = $CMD{coverage};
	$sample1 = $CMD{sample1};
	$sample2 = $CMD{sample2};
	$length = $CMD{length};

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
