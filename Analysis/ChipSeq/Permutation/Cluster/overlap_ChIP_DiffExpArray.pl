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
#  Module: Analysis::ChipSeq::Permutation::Cluster::overlap_ChIP_DiffExpArray.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# Overlap permutation test for
# 1. Array data (differentially expressed genes) and 
# 2. ChIP segments
#
# Written by Stephan Ossowski, 06/20/09
# --------------------------------------------------------------------


use Getopt::Long;
use FindBin;
use DBI;

### Variables
my $db_name;
my $chip_table;
my $diffexp_table;
my $permnum;
my $pref;
my $minpeak;
my $minlen;
my $dfd_early_no_rec_hs_up;
my $dfd_early_no_rec_hs_down;
my $dfd_early_gof_up;
my $dfd_early_gof_down;
my $dfd_late_gof_up;
my $dfd_late_gof_down;


### Get command line options
my %CMD;
GetCom();


### Connect to DB
my $dbh;
&connect_to_db();

### Initialize containers
my %real_result = ();
my %permutation_good = ();
my %permutation_bad  = ();


my %chip_gene = ();

### Load genes near chip
my $q= "SELECT  fbgn
	FROM    $chip_table
	WHERE   max_cov >= $minpeak &&
		segment_length >= $minlen
	ORDER BY begin
";
my $sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {
	$chip_gene{$ref->{fbgn}} = 1;
}
$sth->finish();


### Check overlap with differentially expressed genes
$q = "	SELECT distinct sequence_derived_from
	FROM $diffexp_table
	WHERE   sequence_derived_from like 'FB%' &&
		(
			#dfd_early_no_rec_hs > $dfd_early_no_rec_hs_up || 
			#dfd_early_no_rec_hs between 0.001 and $dfd_early_no_rec_hs_down #||
			#dfd_early_rec_hs > $dfd_early_no_rec_hs_up ||
			#dfd_early_rec_hs between 0.001 and $dfd_early_no_rec_hs_down #||
			dfd_early_gof > $dfd_early_gof_up || 
			dfd_early_gof between 0.001 and $dfd_early_gof_down || 
			dfd_late_gof > $dfd_late_gof_up || 
			dfd_late_gof between 0.001 and $dfd_late_gof_down
		)
";
$sth = $dbh->prepare($q);
$sth->execute();

my $diffExpOverlap = 0;
my $diffExpNoOverlap = 0;
my $counter = 0;
while(my $ref = $sth->fetchrow_hashref()) {
	my $diff_gene = $ref->{sequence_derived_from};
	if(exists $chip_gene{$diff_gene}) {
		$diffExpOverlap++;
	}
	else {
		$diffExpNoOverlap++;
	}
	$counter++;
}


### Check 1000 times random gene
$q = "  SELECT distinct sequence_derived_from
	FROM $diffexp_table
	WHERE sequence_derived_from like 'FB%'
";
$sth = $dbh->prepare($q);
$sth->execute();

my @any_gene = ();

### Store peaks for permutation
while(my $ref = $sth->fetchrow_hashref()) {
	push(@any_gene, $ref->{sequence_derived_from});
}

my $total_genes = @any_gene;

print "$counter\t$total_genes\t$permnum\n";

### Permute
for( my $permute = 1; $permute <= $permnum; $permute++) {
	for(my $i = 1; $i <= $counter; $i++) {
		my $curr_rand = int(rand($total_genes));
		if( exists $chip_gene{$any_gene[$curr_rand]} ) {
			$permutation_good{$permute}++;
		}
		else {
			$permutation_bad{$permute}++;
		}
	}
}



### Print real overlap
open RESULTS, ">$pref.real.txt" or die "Cannot open outfile $pref.real.txt\n";
print RESULTS "$diffExpOverlap\t$diffExpNoOverlap = " . sprintf("%.2f", $diffExpOverlap/($diffExpNoOverlap+$diffExpOverlap)*100) . "%\n";
close RESULTS;


### Print permutation distribution
open PERMOUT, ">$pref.txt" or die "Cannot open outfile $pref.txt\n";
foreach (sort keys %permutation_good) {
	print PERMOUT "$permutation_good{$_}\t$permutation_bad{$_}\t". sprintf("%.2f", $permutation_good{$_}/($permutation_bad{$_}+$permutation_good{$_})*100) ."%\n";
}
close PERMOUT;


### Print permutation options
open OPTIONS, ">$pref.options.txt" or die "Cannot open outfile $pref.options.txt\n";
print OPTIONS 	"DB\t$db_name\nChIP-Table\t$chip_table\nDiffExpArray-Table\t$diffexp_table\n" .
		"Permutations\t$permnum\nMin-Peak\t$minpeak\nMinSegLength\t$minlen\n";
close OPTIONS;


### Call R plot script
my $cmd = "R --slave --vanilla --args " . $diffExpOverlap
		. " $pref.txt $pref.options.txt $pref.hist.pdf < "
		. $FindBin::Bin . "/overlap_Permutation.R";
print STDERR "$cmd\n";
system($cmd);


exit(0);



### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = $db_name;
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}



### Read command line parameters
sub GetCom
{
   my @usage = ("\nUsage: $0

Data:
--db       STRING     Database name
--diffexp  STRING     Array data table: differentially expressed genes 
--chip     STRING     Cluster of TFBS/Enhancer prediction tool <hox_box | cis_analyst>. Several tools can be specified seperated by comma.

Permutation:
--permnum  INT        Number of permutation rounds, default 1000
--pref     STRING     Prefix for all output files, default 'permutation_results'

ChIP constraints:
--minpeak   INT        Minimum peak height for ChIPseq peaks filter, default 6
--minlen    INT        Minimum peak length for ChIPseq peaks filter, default 100

DiffExp Array constraints:
--dfd_early_no_rec_hs_up     DOUBLE    default 1.4
--dfd_early_no_rec_hs_down   DOUBLE    default 0.6
--dfd_early_gof_up           DOUBLE    default 1.4
--dfd_early_gof_down         DOUBLE    default 0.6
--dfd_late_gof_up            DOUBLE    default 1.4
--dfd_late_gof_down          DOUBLE    default 0.6
\n");


        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "db=s", "diffexp=s", "chip=s", "permnum=s", "pref=s", "minpeak=s", "minlen=s",
                        "dfd_early_no_rec_hs_up=s", "dfd_early_gof_up=s", "dfd_late_gof_up=s",
                        "dfd_early_no_rec_hs_down=s", "dfd_early_gof_down=s", "dfd_late_gof_down=s");


        # Database and tables
        die("Please specify database name.\n") unless defined($CMD{db});
        die("Please specify DiffExp array table name.\n") unless defined($CMD{diffexp});
        die("Please specify ChIP table name.\n") unless defined($CMD{chip});
        $db_name           = $CMD{db};
        $diffexp_table     = $CMD{diffexp};
	$chip_table        = $CMD{chip};

	# Permutations
	$permnum = 1000;
	$pref = "permutation_results";
	if(defined $CMD{permnum}) { $permnum = $CMD{permnum}; }
	if(defined $CMD{pref}) { $pref = $CMD{pref}; }

	# ChIP constraints
	$minpeak = 6;
	$minlen = 100;
	if(defined $CMD{minpeak}) { $minpeak = $CMD{minpeak}; }
	if(defined $CMD{minlen})  { $minlen = $CMD{minlen}; }


        # DiffExp Array constraints
        $dfd_early_no_rec_hs_up   = 1.4;
        $dfd_early_no_rec_hs_down = 0.6;
        $dfd_early_gof_up         = 1.6;
        $dfd_early_gof_down       = 0.4;
        $dfd_late_gof_up          = 1.6;
        $dfd_late_gof_down        = 0.4;
        if(defined $CMD{dfd_early_no_rec_hs_up})   { $dfd_early_no_rec_hs_up = $CMD{dfd_early_no_rec_hs_up}; }
        if(defined $CMD{dfd_early_no_rec_hs_down}) { $dfd_early_no_rec_hs_down = $CMD{dfd_early_no_rec_hs_down}; }
        if(defined $CMD{dfd_early_gof_up})         { $dfd_early_gof_up = $CMD{dfd_early_gof_up}; }
        if(defined $CMD{dfd_early_gof_down})       { $dfd_early_gof_down = $CMD{dfd_early_gof_down}; }
        if(defined $CMD{dfd_late_gof_up})          { $dfd_late_gof_up = $CMD{dfd_late_gof_up}; }
        if(defined $CMD{dfd_late_gof_down})        { $dfd_late_gof_down = $CMD{dfd_late_gof_down}; }
}


