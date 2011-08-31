#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $chip_table = shift; # e.g. enriched_dfd_s1_nonrep

my $dbh;
&connect_to_db();

### Query chromosome length
my $q = "SELECT chromosome, chr_alias FROM seq_max ORDER BY chromosome";
my $sth = $dbh->prepare($q);
$sth->execute();

### Parse predictions 
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr       = $ref->{chromosome};
        my $chr_alias = $ref->{chr_alias};
	my %peak = ();

	### Load ChipSeq peaks
        my $q ="SELECT  chromosome, chr_alias, begin, end, segment_length, avg_cov, max_cov
                FROM    $chip_table
                WHERE   chromosome = $chr
                ORDER by chromosome, begin
        ";
        my $sth = $dbh->prepare($q);
        $sth->execute();
	
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$peak{$i} = "$ref->{chromosome}\t$ref->{chr_alias}\t$ref->{begin}\t$ref->{end}\t$ref->{segment_length}\t$ref->{avg_cov}\t$ref->{max_cov}";
		}
	}

	### Load segments (intergenic and intron)
        $q = "	SELECT  fbgn, name, segment_start_position, segment_end_position
                FROM    segment_2_gene
                WHERE   chromosome = '$chr_alias'
                ORDER BY segment_start_position
        ";
        $sth = $dbh->prepare($q);
        $sth->execute();

	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{segment_start_position};
		my $end = $ref->{segment_end_position};
		my %hit = ();

		for(my $pos = $beg; $pos <= $end; $pos++) {
			if(exists $peak{$pos}) { $hit{"$peak{$pos}\t$ref->{fbgn}\t$ref->{name}"} = 1; }
		}
		foreach (sort keys %hit) { 
			print "$_\n";
		}
	}
}

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

