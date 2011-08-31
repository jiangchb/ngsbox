#!/usr/bin/perl

# --------------------------------------------------------------------
# Overlap permutation test for
# 1. ChIPseq peaks, 
# 2. Array data (differentially expressed genes) and 
# 3. TFBS/Enhancer predictions
#
# Written by Stephan Ossowski, 06/20/09
# --------------------------------------------------------------------


use strict;
use warnings;
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
my $minpeak;
my $minlen;
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

### Open outfile for regulated gene candidates
open GENE, ">$pref.gene_candidate.txt" or die "\n\nCannot open outfile.\n\n";

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

	### Load ChIPseq peaks
	my $q= "SELECT distinct begin, end
		FROM    $chip_table
		WHERE   chr = $chr &&
			max_cov >= $minpeak &&
			segment_length >= $minlen
		ORDER by chr, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin} - $range; $i <= $ref->{end} + $range; $i++) {
			$peak{$i} = 1;
		}
	}


        ### Load TFBS prediction
        $q = "	SELECT  begin, end, fbgn
                FROM    cluster_cisA_hoxB_gene
                WHERE   ($tool_selection) &&
			(chromosome = '$chr_alias') &&
			#(averageCons >= 0.3) &&
			(number_BS >= 3) &&
			( (number_species >= $cons) || (number_species IS NULL) ) &&
			fbgn IN(SELECT sequence_derived_from 
				FROM $diffexp_table
				WHERE   dfd_early_no_rec_hs > $dfd_early_no_rec_hs_up || 
					dfd_early_no_rec_hs between 0.001 and $dfd_early_no_rec_hs_down ||
					#dfd_early_rec_hs > $dfd_early_no_rec_hs_up ||
					#dfd_early_rec_hs between 0.001 and $dfd_early_no_rec_hs_down #||
					dfd_early_gof > $dfd_early_gof_up || 
					dfd_early_gof between 0.001 and $dfd_early_gof_down || 
					dfd_late_gof > $dfd_late_gof_up || 
					dfd_late_gof between 0.001 and $dfd_late_gof_down
				)
                ORDER by chromosome, begin
        ";
        $sth = $dbh->prepare($q);
        $sth->execute();

	my %cluster = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};

		### Store peaks for permutation
		$cluster{$beg} = $end;

		### Check real result overlap
		my $hit = 0;
		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $peak{$pos}) { $hit = 1; }
		}
		if($hit == 1) { 
			$real_result{good}++; 
			print GENE "$chr_alias:$beg..$end\t" . $ref->{fbgn} . "\n";
		}
		else          { $real_result{bad}++; }
	}
	$sth->finish();


	### Permute
	for( my $permute = 1; $permute <= $permnum; $permute++) {
		my %result = (good => 0, bad => 0);
		my $curr_rand = int(rand($max_pos)) + 1;
		
		foreach my $beg (keys %cluster) {
			my $end = $cluster{$beg};
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
print RESULTS "$real_result{good}\t$real_result{bad}\t" . sprintf("%.2f", $real_result{good}/($real_result{good}+$real_result{bad})*100) . "%\n";
close RESULTS;


### Print permutation distribution
open PERMOUT, ">$pref.txt" or die "Cannot open outfile $pref.txt\n";
foreach (sort keys %permutation_good) {
	print PERMOUT "$permutation_good{$_}\t$permutation_bad{$_}\t". sprintf("%.2f", $permutation_good{$_}/($permutation_good{$_}+$permutation_bad{$_})*100) ."%\n";
}
close PERMOUT;


### Print permutation options
open OPTIONS, ">$pref.options.txt" or die "Cannot open outfile $pref.options.txt\n";
print OPTIONS 	"DB\t$db_name\nChIP-Table\t$chip_table\nDiffExpArray-Table\t$diffexp_table\n" .
		"Cluster-Table\t$cluster_tool_list\nPermutations\t$permnum\nMin-Peak\t$minpeak\nMinSegLength\t$minlen\nRange\t$range";
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
--diffexp  STRING     Array data table: differentially expressed genes 
--cluster  STRING     Cluster of TFBS/Enhancer prediction tool <hox_box | cis_analyst>. Several tools can be specified seperated by comma.

Permutation:
--permnum  INT        Number of permutation rounds, default 1000
--pref     STRING     Prefix for all output files, default 'permutation_results'

Overlap constraints:
--range     INT        Extend overlap range, default 100

ChIP constraints:
--minpeak   INT        Minimum peak height for ChIPseq peaks filter, default 6
--minlen    INT        Minimum peak length for ChIPseq peaks filter, default 100

Cluster constraints:
--cons      INT        Minimum number of species the cluster is conserved in.

DiffExp Array constraints:
--dfd_early_no_rec_hs_up     DOUBLE    default 1.4
--dfd_early_no_rec_hs_down   DOUBLE    default 0.6
--dfd_early_gof_up           DOUBLE    default 1.4
--dfd_early_gof_down         DOUBLE    default 0.6
--dfd_late_gof_up            DOUBLE    default 1.4
--dfd_late_gof_down          DOUBLE    default 0.6
\n");


        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "db=s", "chip=s", "diffexp=s", "cluster=s", "permnum=s", "pref=s", "minpeak=s", "minlen=s", "range=s", "cons=s",
                        "dfd_early_no_rec_hs_up=s", "dfd_early_gof_up=s", "dfd_late_gof_up=s",
                        "dfd_early_no_rec_hs_down=s", "dfd_early_gof_down=s", "dfd_late_gof_down=s");


        # Database and tables
        die("Please specify database name.\n") unless defined($CMD{db});
        die("Please specify ChIPseq paet table name.\n") unless defined($CMD{chip});
        die("Please specify DiffExp array table name.\n") unless defined($CMD{diffexp});
        die("Please specify bindingsite table name.\n") unless defined($CMD{cluster});
        $db_name           = $CMD{db};
        $chip_table        = $CMD{chip};
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


        # ChIP constraints:
        $minpeak = 6;
        $minlen = 100;
        if(defined $CMD{minpeak}) { $minpeak = $CMD{minpeak}; }
        if(defined $CMD{minlen})  { $minlen = $CMD{minlen}; }

	
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
