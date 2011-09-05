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
#  Module: Analysis::ChipSeq::Permutation::TFBS::overlap_TFBS_DiffExpArray.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# Overlap permutation test for
# 1. Array data (differentially expressed genes) and 
# 2. TFBS/Enhancer predictions
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
my $cluster_tool_list;
my $permnum;
my $pref;
my $range;
my $cons;
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

### Get selected prediction tools
my $tool_selection = "";
my @cluster_tools = split(/,/, $cluster_tool_list);
foreach my $cluster_tool (@cluster_tools) {
	$tool_selection .= "prediction_tool = '$cluster_tool' || ";
}
$tool_selection = substr($tool_selection, 0, length($tool_selection) -3 );



my %cluster_gene = ();
print "HALLO\n\n";
### Load genes near cluster
my $q= "SELECT  fbgn
	FROM    tfbs_cisA_hoxB_gene
	WHERE   ($tool_selection) &&
		(average_conservation > 0.3)
	ORDER BY begin
";
my $sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {
	$cluster_gene{$ref->{fbgn}} = 1;
}
$sth->finish();


### Check overlap with differentially expressed genes
$q = "	SELECT sequence_derived_from
	FROM $diffexp_table
	WHERE   #dfd_early_no_rec_hs > $dfd_early_no_rec_hs_up || 
		#dfd_early_no_rec_hs between 0.001 and $dfd_early_no_rec_hs_down || 
		dfd_early_gof > $dfd_early_gof_up || 
		dfd_early_gof between 0.001 and $dfd_early_gof_down || 
		dfd_late_gof > $dfd_late_gof_up || 
		dfd_late_gof between 0.001 and $dfd_late_gof_down
";
$sth = $dbh->prepare($q);
$sth->execute();

my $diffExpOverlap = 0;
my $diffExpNoOverlap = 0;
my $counter++;
while(my $ref = $sth->fetchrow_hashref()) {
	my $diff_gene = $ref->{sequence_derived_from};
	if(exists $cluster_gene{$diff_gene}) {
		$diffExpOverlap++;
	}
	else {
		$diffExpNoOverlap++;
	}
	$counter++;
}


### Check 1000 times random gene
$q = "  SELECT sequence_derived_from
	FROM $diffexp_table
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
		if( exists $cluster_gene{$any_gene[$curr_rand]} ) {
			$permutation_good{$permute}++;
		}
		else {
			$permutation_bad{$permute}++;
		}
	}
}



### Print real overlap
open RESULTS, ">$pref.real.txt" or die "Cannot open outfile $pref.real.txt\n";
print RESULTS "$diffExpOverlap/$diffExpNoOverlap = " . sprintf("%.2f", $diffExpOverlap/$diffExpNoOverlap*100) . "%\n";
close RESULTS;


### Print permutation distribution
open PERMOUT, ">$pref.txt" or die "Cannot open outfile $pref.txt\n";
foreach (sort keys %permutation_good) {
	print PERMOUT "$permutation_good{$_}\t$permutation_bad{$_}\t". sprintf("%.2f", $permutation_good{$_}/$permutation_bad{$_}*100) ."%\n";
}
close PERMOUT;


### Print permutation options
open OPTIONS, ">$pref.options.txt" or die "Cannot open outfile $pref.options.txt\n";
print OPTIONS 	"DB\t$db_name\nDiffExpArray-Table\t$diffexp_table\n" .
		"Cluster-Table\t$cluster_tool_list\nPermutations\t$permnum\nRange\t$range";
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
--cluster  STRING     Cluster of TFBS/Enhancer prediction tool <hox_box | cis_analyst>. Several tools can be specified seperated by comma.

Permutation:
--permnum  INT        Number of permutation rounds, default 1000
--pref     STRING     Prefix for all output files, default 'permutation_results'

Overlap constraints:
--range     INT        Extend overlap range, default 100

Cluster constraints:
--cons      INT        Minimum number of species the cluster is conserved in.

DiffExp Array constraints:
--dfd_early_no_rec_hs_up     DOUBLE    default 1.5
--dfd_early_no_rec_hs_down   DOUBLE    default 0.5
--dfd_early_gof_up           DOUBLE    default 1.5
--dfd_early_gof_down         DOUBLE    default 0.5
--dfd_late_gof_up            DOUBLE    default 1.5
--dfd_late_gof_down          DOUBLE    default 0.5
\n");


        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "db=s", "diffexp=s", "cluster=s", "permnum=s", "pref=s", "range=s", "cons=s",
                        "dfd_early_no_rec_hs_up=s", "dfd_early_gof_up=s", "dfd_late_gof_up=s",
                        "dfd_early_no_rec_hs_down=s", "dfd_early_gof_down=s", "dfd_late_gof_down=s");


        # Database and tables
        die("Please specify database name.\n") unless defined($CMD{db});
        die("Please specify DiffExp array table name.\n") unless defined($CMD{diffexp});
        die("Please specify bindingsite table name.\n") unless defined($CMD{cluster});
        $db_name           = $CMD{db};
        $diffexp_table     = $CMD{diffexp};
        $cluster_tool_list = $CMD{cluster};


	# Permutations
	$permnum = 1000;
	$pref = "permutation_results";
	if(defined $CMD{permnum}) { $permnum = $CMD{permnum}; }
	if(defined $CMD{pref}) { $pref = $CMD{pref}; }


	# Overlap constrainst
	$range = 100;
	if(defined $CMD{range}) { $range = $CMD{range}; }


	# Cluster constraints:
	$cons = 1;
	if(defined $CMD{cons}) { $cons = $CMD{cons}; }

        # DiffExp Array constraints
        $dfd_early_no_rec_hs_up   = 1.4;
        $dfd_early_no_rec_hs_down = 0.6;
        $dfd_early_gof_up         = 1.4;
        $dfd_early_gof_down       = 0.6;
        $dfd_late_gof_up          = 1.4;
        $dfd_late_gof_down        = 0.6;
        if(defined $CMD{dfd_early_no_rec_hs_up})   { $dfd_early_no_rec_hs_up = $CMD{dfd_early_no_rec_hs_up}; }
        if(defined $CMD{dfd_early_no_rec_hs_down}) { $dfd_early_no_rec_hs_down = $CMD{dfd_early_no_rec_hs_down}; }
        if(defined $CMD{dfd_early_gof_up})         { $dfd_early_gof_up = $CMD{dfd_early_gof_up}; }
        if(defined $CMD{dfd_early_gof_down})       { $dfd_early_gof_down = $CMD{dfd_early_gof_down}; }
        if(defined $CMD{dfd_late_gof_up})          { $dfd_late_gof_up = $CMD{dfd_late_gof_up}; }
        if(defined $CMD{dfd_late_gof_down})        { $dfd_late_gof_down = $CMD{dfd_late_gof_down}; }
}


