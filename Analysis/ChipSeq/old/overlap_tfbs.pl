#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use DBI;

my %CMD = ();

# parameters
my $range 		= 100;
my $max_cov 		= 6;
my $segment_length 	= 100;
my $permutations 	= 1000;
my $conservation	= -1;

# parameters
my @tfbs_tables 	= ();

# output
my %total_good 		= ();
my %total_bad 		= ();
my %result 		= (good => 0, bad => 0);
my %positions		= ();

my $dbh;

### get the parameters
GetCom();
&connect_to_db();

### Query chromosome length
my $q0 = "SELECT chr, chr_alias, max_pos FROM seq_max where chr BETWEEN 3 AND 8 ORDER BY chr";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();

### delete files
unlink('permutations.txt');
unlink('result.txt');
unlink('options.txt');
unlink('positions');

### write positions results
my $filename = "positions_r_".$range."_m_".$max_cov."_s_".$segment_length."_t_".join(",", @tfbs_tables)."_c_".$conservation.".txt";
unlink($filename);
open (MY_POSITIONS_FILE, '>'.$filename);

### For each chromosome
while( my $ref0 = $sth0->fetchrow_hashref() ) {
	my $chr       = $ref0->{'chr'};
	my $chr_alias = $ref0->{chr_alias};
	my $max_pos   = $ref0->{max_pos};
	my %tfbs 	= ();

	### Load TFBS prediction
	### ("binding_site", "cluster", "redfly_tfbs", "redfly_crm")
	foreach my $tfbs_table (@tfbs_tables) {
		my $q = "";
		if($tfbs_table eq "binding_site"){
			$q= "   SELECT 	begin, end
				FROM 	binding_site
				WHERE 	chromosome = '$chr_alias'
				ORDER BY begin
			";
		} elsif($tfbs_table eq "binding_site_hox_box"){
			my $str = "";
			if($conservation >= 0){
				$str = " AND binding_site.patser_score > $conservation ";
			}
			$q= "   SELECT 	binding_site.begin, binding_site.end
				FROM 	binding_site, cluster
				WHERE 	binding_site.chromosome = '$chr_alias' $str AND cluster.old_prediction = 0 AND cluster.id = binding_site.parent_cluster_id
				ORDER BY binding_site.begin
			";
		} elsif($tfbs_table eq "binding_site_cis_analyst"){
			my $str = "";
			if($conservation >= 0){
				$str = " AND binding_site.average_conservation > $conservation ";
			}
			$q= "   SELECT 	binding_site.begin, binding_site.end
				FROM 	binding_site, cluster
				WHERE 	binding_site.chromosome = '$chr_alias' $str AND cluster.old_prediction = 1 AND cluster.id = binding_site.parent_cluster_id
				ORDER BY binding_site.begin
			";
		} elsif ($tfbs_table eq "cluster_hox_box"){
			$q= "   SELECT 	begin, end
				FROM 	cluster_hox_box
				WHERE 	chromosome = '$chr_alias' AND old_prediction = 0
				ORDER BY begin
			";
		} elsif ($tfbs_table eq "cluster_cis_analyst"){
			$q= "   SELECT 	begin, end
				FROM 	cluster_cis_analyst
				WHERE 	chromosome = '$chr_alias' AND old_prediction = 1
				ORDER BY begin
			";
		} elsif ($tfbs_table eq "cluster"){
			$q= "   SELECT 	begin, end
				FROM 	cluster_all
				WHERE 	chromosome = '$chr_alias'
				ORDER BY begin
			";
		} else {
			$q= "   SELECT 	begin, end
				FROM 	$tfbs_table
				WHERE 	chromosome = '$chr_alias'
				ORDER BY begin
			";
		}

		my $sth = $dbh->prepare($q);
		$sth->execute();

		while(my $ref = $sth->fetchrow_hashref()) {
			for(my $i = $ref->{begin} - $range; $i <= $ref->{end} + $range; $i++) {
				$tfbs{$i} = 1;
			}
		}
	}

	### Load ChipSeq peaks
	my $q="	SELECT 	distinct begin, end, name
		FROM	tool_shore_enriched_dfd_50_7_gene
		WHERE 	chr = $chr &&
			max_cov >= $max_cov &&
			segment_length >= $segment_length &&
			fbgn in (select sequence_derived_from from ann_array_dfd_anno_gene_jun08 where dfd_late_gof > 1.5 AND dfd_late_gof BETWEEN 0.001 AND 0.5 AND dfd_early_no_rec_hs > 1.5 || dfd_early_no_rec_hs BETWEEN 0.001 AND 0.5 || dfd_early_gof > 1.5 || dfd_early_gof between 0.001 and 0.5)
		ORDER by chr, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %peak = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};
		$peak{$beg} = $end;
		my $hit = 0;

		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $tfbs{$pos}) { $hit = 1; }
		}
		if($hit == 1) { 
			$result{good}++;
			print MY_POSITIONS_FILE "$chr_alias:$beg..$end\t$ref->{name}\n";
		} else { $result{bad}++; }
	}

	### Compare
	for( my $permute = 0; $permute < $permutations; $permute++) {
		my %result_local = (good => 0, bad => 0);
		my $curr_rand = int(rand($max_pos)) + 1;
		
		foreach my $beg (keys %peak) {
			my $end = $peak{$beg};
			my $rand_beg = (($beg + $curr_rand)%$max_pos) + 1;
			my $rand_end = (($end + $curr_rand)%$max_pos) + 1;
		
			my $hit = 0;

			for(my $pos = $rand_beg; $pos <= $rand_end; $pos++) {
				if(exists $peak{$pos}) { $hit = 1; }
			}
			if($hit == 1) { $result_local{good}++; }
			else          { $result_local{bad}++; }
		}
		$total_good{$permute} += $result_local{good};
		$total_bad{$permute} += $result_local{bad};
	}
}
$sth0->finish();
close (MY_POSITIONS_FILE);

### write permutation results
open (MY_PERMUTATIONS_FILE, '>permutations.txt');
foreach (sort keys %total_good) {
	print MY_PERMUTATIONS_FILE "$total_good{$_}\t$total_bad{$_}\t". sprintf("%.2f", $total_good{$_}/$total_bad{$_}*100) ."%\n";
}
close (MY_PERMUTATIONS_FILE);

### write result
open (MY_RESULT_FILE, '>result.txt');
print MY_RESULT_FILE "$result{good}\t$result{bad}\t" . sprintf("%.2f", $result{good}/$result{bad}*100) . "%\n";
close (MY_RESULT_FILE);

### write options
open (MY_OPTIONS_FILE, '>options.txt');
my $string = join(",", @tfbs_tables);
print MY_OPTIONS_FILE "range\t$range\nmax_cov\t$max_cov\nsegment_length\t$segment_length\npermutations\t$permutations\ntfbs_tables\t$string\nconservation\t$conservation\n";
close (MY_OPTIONS_FILE);

### call R
system("R --slave --vanilla < overlap_tfbs.R");

exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "chip_seq_dmel";
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

sub GetCom {

  my @usage = ("\nUsage: $0

required:
--tfbs_tables\t\tcomma separeted: binding_site,cluster, .... (add _hox_box or _cis_analyst)
\n
optional:
--max_cov\t\tstandard = $max_cov
--segment_length\tstandard = $segment_length
--range\t\t\tstandard = $range
--permutations\t\tstandard = $permutations
--conservation\t\tdepends on your table _hox_box knows patser_score(0-16) and _cis_analyst knows average_conservation (0-1); otherwise omitted.
\n");

	die(@usage) if (@ARGV == 0);
        
        GetOptions(\%CMD, "tfbs_tables=s", "max_cov=s", "segment_length=s", "range=s", "permutations=s", "conservation=s");

	die(@usage) if (!defined($CMD{tfbs_tables}));

	my $string 	= $CMD{tfbs_tables} if(defined($CMD{tfbs_tables}));
	@tfbs_tables	= split(",", $string);
        $max_cov        = $CMD{max_cov} if(defined($CMD{max_cov}));
        $segment_length = $CMD{segment_length} if(defined$CMD{segment_length});
        $range 		= $CMD{range} if(defined($CMD{range}));
        $permutations 	= $CMD{permutations} if(defined($CMD{permutations}));
	$conservation	= $CMD{conservation} if(defined($CMD{conservation}));

        return(0);
}

