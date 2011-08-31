#!/usr/bin/perl

# --------------------------------------------------------------------
# Overlap and Permutation test for 
# 1. ChIPseq peaks, 
# 2. Array data (differentially expressed genes) and 
# 3. TFBS cluster/Enhancer predictions
#
# Written by Stephan Ossowski, 06/20/09
# --------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use DBI;


### Variables
my $db_name;
my $chip_table;
my $diffexp_table;
my $bs_table;
my $minpeak;
my $minlen;
my $range;
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
my @tfbs_tables = split(/,/, $bs_table);
my %result = ();


### Query chromosome length
my $q = "SELECT chr, chr_alias FROM seq_max where chr between 3 and 8 ORDER BY chr";
my $sth = $dbh->prepare($q);
$sth->execute();

### Parse TFBS predictions 
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr       = $ref->{'chr'};
        my $chr_alias = $ref->{chr_alias};
	my %tfbs = ();


	### Load TFBS prediction
	foreach my $tfbs_table (@tfbs_tables) {
		my $q= "SELECT 	begin, end
			FROM 	$tfbs_table
			WHERE 	chromosome = '$chr_alias'
			ORDER BY begin
		";
		my $sth = $dbh->prepare($q);
		$sth->execute();
	
		while(my $ref = $sth->fetchrow_hashref()) {
			for(my $i = $ref->{begin} - $range; $i <= $ref->{end} + $range; $i++) {
				$tfbs{$i} = 1;
			}
		}
	}


	### Load ChipSeq peaks
	my $q="	SELECT 	distinct begin, end
		FROM	$chip_table
		WHERE 	chr = $chr &&
			max_cov >= $minpeak &&
			segment_length >= $minlen &&
			fbgn IN(SELECT sequence_derived_from 
				FROM $diffexp_table
				WHERE 	dfd_early_no_rec_hs > $dfd_early_no_rec_hs_up || 
					dfd_early_no_rec_hs between 0.001 and $dfd_early_no_rec_hs_down || 
					dfd_early_gof > $dfd_early_gof_up || 
					dfd_early_gof between 0.001 and $dfd_early_gof_down || 
					dfd_late_gof > $dfd_late_gof_up || 
					dfd_late_gof between 0.001 and $dfd_late_gof_down)
		ORDER by chr, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};
		my $hit = 0;

		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $tfbs{$pos}) { $hit = 1; }
		}
		if($hit == 1) { $result{good}++; print "$chr_alias:$beg..$end\n";}
		else          { $result{bad}++; }
	}
}

print "$result{good}/$result{bad} = " . sprintf("%.2f", $result{good}/$result{bad}*100) . "%\n";

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
--bs       STRING     Bindingsite/Enhancer table. Several table can be specified at one (comma seperated).

Overlap constraints
--range     INT        Extend overlap range, default 100

ChIP constraints:
--minpeak   INT        Minimum peak height for ChIPseq peaks filter, default 6
--minlen    INT        Minimum peak length for ChIPseq peaks filter, default 100

DiffExp Array constraints
--dfd_early_no_rec_hs_up     DOUBLE    default 1.5
--dfd_early_no_rec_hs_down   DOUBLE    default 0.5
--dfd_early_gof_up           DOUBLE    default 1.5
--dfd_early_gof_down         DOUBLE    default 0.5
--dfd_late_gof_up            DOUBLE    default 1.5
--dfd_late_gof_down          DOUBLE    default 0.5
\n");


        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "db=s", "chip=s", "diffexp=s", "bs=s", "minpeak=s", "minlen=s", "range=s",
                        "dfd_early_no_rec_hs_up=s", "dfd_early_gof_up=s", "dfd_late_gof_up=s",
                        "dfd_early_no_rec_hs_down=s", "dfd_early_gof_down=s", "dfd_late_gof_down=s");


        # Database and tables
        die("Please specify database name.\n") unless defined($CMD{db});
        die("Please specify ChIPseq paet table name.\n") unless defined($CMD{chip});
        die("Please specify DiffExp array table name.\n") unless defined($CMD{diffexp});
        die("Please specify bindingsite table name.\n") unless defined($CMD{bs});
        $db_name       = $CMD{db};
        $chip_table    = $CMD{chip};
        $diffexp_table = $CMD{diffexp};
        $bs_table      = $CMD{bs};


        # ChIP constraints:
        $minpeak = 6;
        $minlen = 100;
        if(defined $CMD{minpeak}) { $minpeak = $CMD{minpeak}; }
        if(defined $CMD{minlen})  { $minlen = $CMD{minlen}; }


	# Overlap constrainst
	$range = 100;
	if(defined $CMD{range}) { $range = $CMD{range}; }


        # DiffExp Array constraints
        $dfd_early_no_rec_hs_up   = 1.5;
        $dfd_early_no_rec_hs_down = 0.5;
        $dfd_early_gof_up         = 1.5;
        $dfd_early_gof_down       = 0.5;
        $dfd_late_gof_up          = 1.5;
        $dfd_late_gof_down        = 0.5;
        if(defined $CMD{dfd_early_no_rec_hs_up})   { $dfd_early_no_rec_hs_up = $CMD{dfd_early_no_rec_hs_up}; }
        if(defined $CMD{dfd_early_no_rec_hs_down}) { $dfd_early_no_rec_hs_down = $CMD{dfd_early_no_rec_hs_down}; }
        if(defined $CMD{dfd_early_gof_up})         { $dfd_early_gof_up = $CMD{dfd_early_gof_up}; }
        if(defined $CMD{dfd_early_gof_down})       { $dfd_early_gof_down = $CMD{dfd_early_gof_down}; }
        if(defined $CMD{dfd_late_gof_up})          { $dfd_late_gof_up = $CMD{dfd_late_gof_up}; }
        if(defined $CMD{dfd_late_gof_down})        { $dfd_late_gof_down = $CMD{dfd_late_gof_down}; }
}

