#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $file = shift;
my $callability_check = shift;

if(! defined $callability_check) { $callability_check = 0; }

my $dbh;
&connect_to_db();

open FILE, $file or die "\n\nCannot open infile\n\n";
while( <FILE> ) {
	my @a = split " ", $_;
	my $strain         = $a[0];
	my $chr            = $a[1];
	my $beg_variant    = $a[2];
	my $end_variant    = $a[3];
	my $variant        = $a[4];

	my $beg_sanger  = $beg_variant - 500;
	my $end_sanger  = $end_variant + 500;


	### Get reference sequence
	my $q= "SELECT	position, base
		FROM 	seq_ref
		WHERE	chromosome = $chr && position BETWEEN $beg_sanger AND $end_sanger
		ORDER BY chromosome, position\n";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my %cons_seq = ();
	my %good_ref = ();
	while( my $ref = $sth->fetchrow_hashref() ) {
		$cons_seq{$ref->{position}} = $ref->{base};
		$good_ref{$ref->{position}} = 0;
	}


	### Check reference call for all lines
	foreach my $table ("poly_reference_30_119", "poly_reference_30_29", "poly_reference_30_49", "poly_reference_30_59", "poly_reference_30_69") {

		### Get consensus sequence of surrounding region
		my $q2="SELECT 	position
			FROM 	$table
			WHERE	chromosome = $chr &&
				position BETWEEN $beg_sanger AND $end_sanger
			ORDER by chromosome, position
		";
		my $sth2 = $dbh->prepare($q2);
		$sth2->execute();

		while( my $ref2 = $sth2->fetchrow_hashref() ) {
			$good_ref{$ref2->{position}}++;
		}
	}


	### Print header
	print ">$chr-$beg_sanger-$end_sanger(VAR=$strain-$beg_variant-$end_variant-$variant):TARGET=450,101\n";

	### Print sequence
	for(my $i = $beg_sanger; $i <= $end_sanger; $i++) {
		
		if( ($callability_check == 0) || (exists$good_ref{$i} && $good_ref{$i} == 5) ) {
			print $cons_seq{$i};
		}
		else { 
			if ($i >= $beg_variant && $i <= $end_variant) {
				print $cons_seq{$i};
			}
			else {
				print "N"; 
			}
		}
	}
	print "\n\n";
}

close FILE;

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ma";
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

