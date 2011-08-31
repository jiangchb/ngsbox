#!/usr/bin/perl
####################################################################################
# Author 	Korbinian Schneeberger, Stephan Ossowski
# Date 		06/15/08
# Version	0.2
# Function	Finds specific upregulated segments in sRNA expression data
####################################################################################

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use DBI;

my $database = "srna";
my $up;
my @down;
my $species;
my $fold;
my $prp;
my $outputfilename;
my $min_cov;
my $avg_hits;

my $dbh;
&connect_to_db();

my %CMD = ();
GetCom();

######### Check number of segments ################
my $seg_num = 0;
my $seg_count = 0;
my $q= "SELECT COUNT(*) AS num FROM segment_$up 
	WHERE species = $species AND avg_coverage >= $min_cov AND avg_avg_hits <= $avg_hits";
my $sth = $dbh->prepare($q);
$sth->execute();
while (my $res = $sth->fetchrow_hashref()) { $seg_num = $res->{num}; }

open FILE, "> $outputfilename" or die "Cannot open $outputfilename\n";


######## Check the segments #######################
$q = "	SELECT chromosome, begin, end, avg_coverage FROM segment_$up 
	WHERE species = $species AND avg_coverage >= $min_cov AND avg_avg_hits <= $avg_hits
	ORDER BY chromosome, begin";
$sth = $dbh->prepare($q);
$sth->execute();
while (my $res = $sth->fetchrow_hashref()) {

	### Show progress
	$seg_count++;
	print $seg_count, "/", $seg_num, "\n" if $seg_count%100 == 0;

	### Initialize variables
	my $chr = $res->{chromosome};	
	my $begin = $res->{begin};
	my $end = $res->{end};
	my $cov = $res->{avg_coverage};
	my $found = 0;
	my $prped_out = 0;

	### Check polymorphism table if defined
	if ($prp ne "NA") {
		my $q2="SELECT 	* FROM prp_$prp 
			WHERE 	(begin BETWEEN $begin AND $end) OR 
				(end BETWEEN $begin AND $end) OR 
				(begin <= $begin and end >= $end)";
		my $sth2 = $dbh->prepare($q2); 
		$sth2->execute();
		while (my $res = $sth2->fetchrow_hashref()) { 
			$prped_out = 1;
		}
	}


	### Check fold change
	if ($prp eq "NA" or $prped_out == 0) {
		DOWNS: foreach my $down_sample (@down) {

			my $q1="SELECT 	chromosome, begin, end FROM segment_$down_sample
				WHERE 	(species = $species) AND
					($cov/avg_coverage < $fold) AND 
					(chromosome = $chr) AND 
					(
					  (begin between $begin AND $end) OR 
					  (end between $begin and $end) OR 
					  (begin <= $begin and end >= $end)
					)
				LIMIT 1";
			my $sth1 = $dbh->prepare($q1); $sth1->execute();
			while (my $res = $sth1->fetchrow_hashref()) {
				$found = 1;
				last DOWNS;
			}
		}

		if ($found == 0) {
			print FILE "Chr",$chr, ":", $begin, "..", $end, "\t", $cov, "\n";
		}
	}
}



### Parse input parameters
sub GetCom {

  my @usage = ("\nUsage: $0 

required:
--up\t\tSet one specific sample id for upregulation
--down\t\tSet one or more sample ids for downregulation (comma separated).
--species\t\tsRNA species (length)
--fold\t\tMinimal fold change between samples.
--min_cov\t\tMinimum max. coverage of a segment in the upregulated sample.

optional:
--prp\t\tSet polymorphism table used to mask regions from comparison
--avg_hits\t\tAverage number of hits the reads in this segment are allowed to have. (Default = 1.9)
-show\t\tSet this flag to see all available sample ids. (Program will exit afterwards.)

Examples:
perl compare_segments.pl -show
perl compare_segments.pl --up=005_0030_3_21 --down=003_0029_4_21,001_0029_2_21 --fold=5 --min_cov=10
perl compare_segments.pl --up=005_0030_3_21 --down=003_0029_4_21,001_0029_2_21 --fold=5 --min_cov=5 -prp=cvi
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "up=s", "down=s", "species=s", "fold=s", "prp=s", "show", "min_cov=s", "avg_hits=s");


	### Show available samples/libraries
	if (defined($CMD{show})) {
		my $q = "SELECT table_name,sample,flowcell,lane,reads_24,assay,organ,biorep,techrep,temperature FROM library_translation";
		my $sth = $dbh->prepare($q); 
		$sth->execute();
		print "Name\tSample\tFlowcell\tLane\tAssay\tOrgan\tTemperature\tBio-Rep\tTech-Rep\n";
		while (my $res = $sth->fetchrow_hashref()) {
			$res->{table_name} =~ s/consensus_//g; 
		        print "$res->{table_name}\t$res->{sample}\t$res->{flowcell}\t$res->{lane}\t$res->{assay}\t$res->{organ}\t$res->{biorep}\t$res->{techrep}\t$res->{temperature}\n";
		}
		print "\n";
		exit(0);	
	}


	### Check input
	die("Please specify sample with upregulation\n") unless $CMD{up};
	die("Please specify samples with downregulations\n") unless $CMD{down};
	die("Please specify sRNA species (length)\n") unless $CMD{species};
	die("Please specify required fold change\n") unless $CMD{fold};
	die("Please specify min coverage\n") unless $CMD{min_cov};

	$up = $CMD{up};
	my $down_string = $CMD{down};
	$species = $CMD{species};
        $fold = $CMD{fold};
	$min_cov = $CMD{min_cov};

	if (defined($CMD{prp})) {
		$prp = $CMD{prp};
	}
	else {
		$prp = "NA";
	}
	
	if (defined($CMD{avg_hits})) {
		$avg_hits = $CMD{avg_hits};
	}
	else {
		$avg_hits = 1.9;
	}

	@down = split ",", $down_string;
	$outputfilename = $up."-vs";
	foreach my $d (@down) {
		$outputfilename .= "-".$d;
	}
	$outputfilename .= "-mincov".$min_cov."-fold".$fold."-prp_".$prp.".txt";
	
	return(0);
}



### Connect to MySQL database
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
