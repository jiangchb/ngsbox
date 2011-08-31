#!/usr/bin/perl

use strict;
use warnings;
use DBI;

my $ecotype    = shift;

### Connect to database
my $dbh;
&connect_to_db();


for(my $chr = 1; $chr <= 5; $chr++) {

	### Load Unsequenced regions
	my $q ="SELECT 	begin, end
		FROM 	unsequenced
		WHERE 	sample = '$ecotype' &&
			chromosome = $chr
		ORDER BY chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %unseq = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq{$i} = 1;
		}
	}

	### Load core unsequenced regions
	$q = "	SELECT  begin, end
		FROM    unsequenced_core
		WHERE   sample = '$ecotype' &&
			chromosome = $chr
		ORDER BY chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %unseq_core = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq_core{$i} = 1;
		}
	}

	
	### Load core nonrep unsequenced regions
	$q = "  SELECT  begin, end
		FROM    unsequenced_cn
		WHERE   sample = '$ecotype' &&
			chromosome = $chr
		ORDER BY chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %unseq_cn = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$unseq_cn{$i} = 1;
		}
	}


	### Load SV inversion
	$q = "	SELECT 	cluster1_region1_begin, cluster2_region2_end, gap1, gap2
		FROM	sv_inversion
		WHERE 	sample = '$ecotype' &&
			chromosome = $chr #&&
			#(gap1 > 0 || gap2 > 0)
		ORDER by chromosome, cluster1_region1_begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %sv = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my @features = ( $ref->{cluster2_region2_end}, $ref->{gap1}, $ref->{gap2} );
		$sv{$ref->{cluster1_region1_begin}} = \@features;
	}


	### Load orphan reads
	$q = "  SELECT  begin, end, pattern_support
		FROM    sv_missing_ends
		WHERE   sample = '$ecotype' &&
			chromosome = $chr
		ORDER by chromosome, begin
	";
	$sth = $dbh->prepare($q);
	$sth->execute();

	my %orphans = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
			$orphans{$i} = $ref->{pattern_support};
		}
	}


	### Compare
	my $total_support_cn = 0;
	my $total_count   = 0;
		
	foreach my $sv_beg (sort {$a<=>$b} keys %sv) {
		my $sv_end = $sv{$sv_beg}[0];
		$total_count++;

		# Unseq support for Insertion
		my $support_core = 0;
		my $support_cn = 0;
		for(my $pos = $sv_beg - 100; $pos <= $sv_end + 100; $pos++) {
			if(exists $unseq_core{$pos}) {
				$support_core = 1;
			}
			if(exists $unseq_cn{$pos}) {
				$support_cn = 1;
			}
		}

		# Orphan support for insertions
		my $orphan_support = 0;
		for(my $pos = $sv_beg; $pos <= $sv_end; $pos++) {
			if(exists $orphans{$pos}) {
				$orphan_support = $orphans{$pos};
				last;
			}
		}

		if(($support_cn == 1) || ($orphan_support > 0) ) { $total_support_cn++; }

		print "$ecotype\t$chr\t$sv_beg\t" . $sv{$sv_beg}[0] . "\t" . $sv{$sv_beg}[1] . "\t" . $sv{$sv_beg}[2] . "\t$support_core\t$support_cn\t$orphan_support\n";
	}

	if($total_count != 0) {
		print STDERR "$chr\t$total_support_cn\t$total_count\t" . $total_support_cn/$total_count . "\n";
	}
}
exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ath_pe";
        my $driver = "mysql";
        my $host = "ume.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

