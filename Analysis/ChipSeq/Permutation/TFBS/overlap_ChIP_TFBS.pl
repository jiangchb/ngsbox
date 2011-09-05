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
#  Module: Analysis::ChipSeq::Permutation::TFBS::overlap_ChIP_TFBS.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# Overlap permutation test for
# 1. ChIPseq peaks, 
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
my $tfbs_tool_list;
my $permnum;
my $pref;
my $minpeak;
my $minlen;
my $range;
my $cons;
my $pats;
my $rep;


### Get command line options
my %CMD;
GetCom();


### Connect to DB
my $dbh;
&connect_to_db();


### Open outfile for regulated gene candidates
open GENE, ">$pref.gene_candidate.txt" or die "\n\nCannot open outfile.\n\n";


### Initialize containers
my %real_result = ();
my %permutation_good = ();
my %permutation_bad  = ();


### Get selected prediction tools
my $tool_selection = "";
my @tfbs_tools = split(/,/, $tfbs_tool_list);
foreach my $tfbs_tool (@tfbs_tools) {
	$tool_selection .= "prediction_tool = '$tfbs_tool' || ";
}
$tool_selection = substr($tool_selection, 0, length($tool_selection) -3 );


### Query chromosome length
my $q0 = "SELECT chr, chr_alias, max_pos FROM ann_chromosome_length where chr BETWEEN 3 AND 8 ORDER BY chr";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();


### For each chromosome
while( my $ref0 = $sth0->fetchrow_hashref() ) {
	my $chr       = $ref0->{'chr'};
	my $chr_alias = $ref0->{chr_alias};
	my $max_pos   = $ref0->{max_pos};
	my %peak = ();

	### Load ChIP peaks
	my $q=" SELECT distinct begin, end
		FROM    $chip_table
		WHERE   chr = $chr &&
			max_cov >= $minpeak &&
			segment_length >= $minlen
		ORDER BY begin";
	my $sth = $dbh->prepare($q);
	$sth->execute();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin} - $range; $i <= $ref->{end} + $range; $i++) {
			$peak{$i} = 1;
		}
	}
	$sth->finish();


        ### Load predicted transcription factor binding sites
        $q = " SELECT distinct begin, end
                FROM    tfbs_cisA_hoxB_gene
                WHERE   ($tool_selection) &&
			chromosome = '$chr_alias' &&
			type = 'dfd' &&
			average_conservation >= $cons &&
			repetitiveness <= $rep &&
			( 
			  (prediction_tool = 'hox_box' && patser_score >= $pats) ||
			  (prediction_tool = 'cis_analyst')
			)
                ORDER by chromosome, begin
        ";
        $sth = $dbh->prepare($q);
        $sth->execute();

	my %tfbs = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};

		### Store TFBS for permutation
		$tfbs{$beg} = $end;

		### Check real result overlap
		my $hit = 0;
		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $peak{$pos}) { 
				$hit = 1; 
			}
		}
		if($hit == 1) { 
			$real_result{good}++; 
			print GENE "$chr_alias:$beg..$end\n";
		}
		else { 
			$real_result{bad}++;
		}
	}
	$sth->finish();


	### Permute
	for( my $permute = 1; $permute <= $permnum; $permute++) {
		my %result = (good => 0, bad => 0);
		my $curr_rand = int(rand($max_pos)) + 1;
		
		foreach my $beg (keys %tfbs) {
			my $end = $tfbs{$beg};
			my $rand_beg = (($beg + $curr_rand)%$max_pos) + 1;
			my $rand_end = (($end + $curr_rand)%$max_pos) + 1;
			my $hit = 0;

			for(my $pos = $rand_beg; $pos <= $rand_end; $pos++) {
				if(exists $peak{$pos}) { $hit = 1; }
			}
			if($hit == 1) { $result{good}++; }
			else          { $result{bad}++; }
		}
		$permutation_good{$permute} += $result{good};
		$permutation_bad{$permute} += $result{bad};
	}
}
$sth0->finish();
close GENE;



### Print real overlap
open RESULTS, ">$pref.real.txt" or die "Cannot open outfile $pref.real.txt\n";
print RESULTS "$real_result{good}\t$real_result{bad}\t" . sprintf("%.2f", $real_result{good}/($real_result{bad}+$real_result{good})*100) . "%\n";
close RESULTS;


### Print permutation distribution
open PERMOUT, ">$pref.txt" or die "Cannot open outfile $pref.txt\n";
foreach (sort keys %permutation_good) {
	print PERMOUT "$permutation_good{$_}\t$permutation_bad{$_}\t". sprintf("%.2f", $permutation_good{$_}/($permutation_bad{$_}+$permutation_good{$_})*100) ."%\n";
}
close PERMOUT;


### Print permutation options
open OPTIONS, ">$pref.options.txt" or die "Cannot open outfile $pref.options.txt\n";
print OPTIONS 	"DB\t$db_name\nChIP-Table\t$chip_table\nTFBS-Table\t$tfbs_tool_list\n" .
		"Permutations\t$permnum\nMin-Peak\t$minpeak\nMinSegLength\t" .
		"$minlen\nRange\t$range\tConservation\t$cons\nPatser-Score\t$pats\nRepetitiveness\t$rep\n";
close OPTIONS;


### Call R plot script
my $cmd = "R --slave --vanilla --args " . $real_result{good}
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
--chip     STRING     ChIPseq peak table
--tfbs     STRING     Bindingsite/Enhancer prediction tool <hox_box | cis_analyst>. Several tools can be specified seperated by comma.

Permutation:
--permnum  INT        Number of permutation rounds, default 1000
--pref     STRING     Prefix for all output files, default 'permutation_results'

Overlap constraints:
--range     INT        Extend overlap range, default 100

ChIP constraints:
--minpeak   INT        Minimum peak height for ChIPseq peaks filter, default 6
--minlen    INT        Minimum peak length for ChIPseq peaks filter, default 100

TFBS constraints:
--cons      DOUBLE     Minimum average conservation of bindingsites, default 0.0
--pats      DOUBLE     Minimum patser score (only HoxBox), default 0.0
--rep       DOUBLE     Maximum repetitiveness ( 0.0 to 1.0, 0.0 = absolutely not repetitive ), default 1.0

\n");


        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "db=s", "chip=s", "tfbs=s", "permnum=s", "pref=s",
			"minpeak=s", "minlen=s", "range=s", "cons=s", "pats=s", "rep=s", "range=s");


        # Database and tables
        die("Please specify database name.\n") unless defined($CMD{db});
        die("Please specify ChIPseq paet table name.\n") unless defined($CMD{chip});
        die("Please specify bindingsite table name.\n") unless defined($CMD{tfbs});
        $db_name        = $CMD{db};
        $chip_table     = $CMD{chip};
        $tfbs_tool_list = $CMD{tfbs};


	# Permutations
	$permnum = 1000;
	$pref = "permutation_results";
	if(defined $CMD{permnum}) { $permnum = $CMD{permnum}; }
	if(defined $CMD{pref}) { $pref = $CMD{pref}; }


	# Overlap constrainst
	$range = 100;
	if(defined $CMD{range}) { $range = $CMD{range}; }


        # ChIP constraints:
        $minpeak = 6;
        $minlen = 100;
        if(defined $CMD{minpeak}) { $minpeak = $CMD{minpeak}; }
        if(defined $CMD{minlen})  { $minlen = $CMD{minlen}; }


	# TFBS constraints:
	$cons = 0.0;
	$pats = 0.0;
	$rep  = 1.0;
	if(defined $CMD{cons}) { $cons = $CMD{cons}; }
	if(defined $CMD{pats}) { $pats = $CMD{pats}; }
	if(defined $CMD{rep})  { $rep  = $CMD{rep};  }
}
