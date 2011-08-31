#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $dbh;
&connect_to_db();

### Get SNPs
my $q ="SELECT chromosome, position, (position+length-1) as end, seq, ecotype 
	FROM poly_small_insertion_singletons
	ORDER by chromosome, position";
my $sth = $dbh->prepare($q);
$sth->execute();

while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr    = $ref->{chromosome};
	my $pos    = $ref->{position};
	my $pos2   = $ref->{end};
	my $start  = $pos - 500;
	my $end    = $pos2 + 500;
	my $seq    = $ref->{seq};
	my $ma     = $ref->{ecotype};

	# Get reference sequence
	my $q2="SELECT	position, base 
		FROM 	seq_ref
		WHERE	chromosome = $chr && position BETWEEN $start AND $end
		ORDER BY chromosome, position\n";
	my $sth2 = $dbh->prepare($q2);
	$sth2->execute();

	my %cons_seq = ();
	my %good_ref = ();
	while( my $ref2 = $sth2->fetchrow_hashref() ) {
		$cons_seq{$ref2->{position}} = $ref2->{base};
		$good_ref{$ref2->{position}} = 0;
	}

	### Check reference call for all lines
	foreach my $table ("poly_reference_30_119", "poly_reference_30_29", "poly_reference_30_49",
				"poly_reference_30_59", "poly_reference_30_69") {

		### Get consensus sequence of surrounding region
		my $q3="SELECT 	position
			FROM 	$table
			WHERE	chromosome = $chr &&
				position BETWEEN $start AND $end
			ORDER by chromosome, position
		";
		my $sth3 = $dbh->prepare($q3);
		$sth3->execute();

		while( my $ref3 = $sth3->fetchrow_hashref() ) {
			$good_ref{$ref3->{position}}++;
		}
	}

	# Print sequence for primer design
	print ">$chr-$start-$end(Indel=MA_$ma-$pos-$pos2-$seq):TARGET=450,101\n";

	for(my $i = $start; $i <= $end; $i++) {
		if( (exists$good_ref{$i}) && ($good_ref{$i} == 5) ) {
			print $cons_seq{$i};
		}
		else { 
			#print "N"; 
			print $cons_seq{$i};
		}
	}
	print "\n\n";
}

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ma";
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

