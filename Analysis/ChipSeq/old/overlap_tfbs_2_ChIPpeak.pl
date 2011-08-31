#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $peak           = shift;
my $enriched_table = shift;
my $tfbs_in        = shift;

my @tfbs_tables = split(/,/, $tfbs_in);

my %result = ();

my $dbh;
&connect_to_db();

### Query chromosome length
my $q = "SELECT chromosome, chr_alias FROM seq_max where chromosome between 3 and 8 ORDER BY chromosome";
my $sth = $dbh->prepare($q);
$sth->execute();

### Parse TFBS predictions 
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr       = $ref->{chromosome};
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
			for(my $i = $ref->{begin} - 100; $i <= $ref->{end} + 100; $i++) {
				$tfbs{$i} = 1;
			}
		}
	}

	### Load ChipSeq peaks
	my $q="	SELECT 	distinct begin, end
		FROM	$enriched_table
		WHERE 	chromosome = $chr &&
			max_cov >= $peak &&
			segment_length >= 100 &&
			fbgn in (select sequence_derived_from from array_dfd_anno_gene_jun08 where dfd_early_no_rec_hs > 1.5 || dfd_early_no_rec_hs between 0.001 and 0.5 || dfd_early_gof > 1.5 || dfd_early_gof between 0.001 and 0.5 || dfd_late_gof > 1.5 || dfd_late_gof between 0.001 and 0.5)
		ORDER by chromosome, begin
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

#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "chip_seq_dmel";
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

