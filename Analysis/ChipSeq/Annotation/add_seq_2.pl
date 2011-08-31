#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $dbh;
&connect_to_db();

my $dbh2;
&connect_to_db2();

### Load binding site table
my $q ="SELECT  grh_chr,begin_grh,end_grh,strand_grh,seq_grh,begin_500,end_500,
		chr_fosD,begin_fosD,end_fosD,strand_fosD,seq_fosD,
		chr_ETS,begin_ETS,end_ETS,strand_ETS,seq_ETS,
		name,gene_start_position,gene_end_position,gene_strand,location_type,gene_id,segment_id
	FROM	fosD2ETS2grh2segment2gene2cons
	ORDER by grh_chr, begin_grh
	";
	my $sth = $dbh2->prepare($q);
	$sth->execute();
	
### Add average conservation and repeat probability for each TFBS
while(my $ref = $sth->fetchrow_hashref()) {
	my $chr1 = $ref->{grh_chr};
	my $beg1 = $ref->{begin_grh};
	my $end1 = $ref->{end_grh};
	my $length1 = $end1 - $beg1 + 1;

	### Load conservation, repeat status and sequence for first binding site
        my $q1="SELECT  conservation, isRepeat, nucleotide
                FROM    ann_sequence
                WHERE   chromosome = '$chr1' AND
			position BETWEEN $beg1 AND $end1
                ORDER BY chromosome, position
        ";
        my $sth1 = $dbh->prepare($q1);
        $sth1->execute();

	my $cons_sum1 = 0;
	my $rep_sum1 = 0;
	my $seq1 = "";
	while(my $ref1 = $sth1->fetchrow_hashref()) {
		if(defined $ref1->{conservation}) {
			$cons_sum1 += $ref1->{conservation};
		}
		else { 
			$cons_sum1 += 0;
		}
		$rep_sum1 += $ref1->{isRepeat};
		$seq1 .= $ref1->{nucleotide};
	}
	my $avg_cons1 = $cons_sum1 / $length1;
	my $avg_rep1  = $rep_sum1 / $length1;


	### Load conservation, repeat status and sequence for second binding site
        my $chr2 = $ref->{chr_fosD};
        my $beg2 = $ref->{begin_fosD};
        my $end2 = $ref->{end_fosD};
        my $length2 = $end2 - $beg2 + 1;

        my $q2="SELECT  conservation, isRepeat, nucleotide
                FROM    ann_sequence
                WHERE   chromosome = '$chr2' AND
                        position BETWEEN $beg2 AND $end2
                ORDER BY chromosome, position
        ";
        my $sth2 = $dbh->prepare($q2);
        $sth2->execute();

        my $cons_sum2 = 0;
        my $rep_sum2 = 0;
        my $seq2 = "";
        while(my $ref2 = $sth2->fetchrow_hashref()) {
                if(defined $ref2->{conservation}) {
                        $cons_sum2 += $ref2->{conservation};
                }
                else { 
                        $cons_sum2 += 0;
                }
                $rep_sum2 += $ref2->{isRepeat};
                $seq2 .= $ref2->{nucleotide};
        }
        my $avg_cons2 = $cons_sum2 / $length2;
        my $avg_rep2  = $rep_sum2 / $length2;


	### Load conservation, repeat status and sequence for third binding site
        my $chr3 = $ref->{'chr_ETS'};
        my $beg3 = $ref->{begin_ETS};
        my $end3 = $ref->{end_ETS};
        my $length3 = $end3 - $beg3 + 1;

        my $q3="SELECT  conservation, isRepeat, nucleotide
                FROM    ann_sequence
                WHERE   chromosome = '$chr3' AND
                        position BETWEEN $beg3 AND $end3
                ORDER BY chromosome, position
        ";
        my $sth3 = $dbh->prepare($q3);
        $sth3->execute();

        my $cons_sum3 = 0;
        my $rep_sum3 = 0;
        my $seq3 = "";
        while(my $ref3 = $sth3->fetchrow_hashref()) {
                if(defined $ref3->{conservation}) {
                        $cons_sum3 += $ref3->{conservation};
                }
                else {
                        $cons_sum3 += 0;
                }
                $rep_sum3 += $ref3->{isRepeat};
                $seq3 .= $ref3->{nucleotide};
        }
        my $avg_cons3 = $cons_sum3 / $length3;
        my $avg_rep3  = $rep_sum3 / $length3;


	print 	"$chr1\t$beg1\t$end1\t" . $ref->{strand_grh} ."\t". $ref->{seq_grh} ."\t". sprintf("%.3f", $avg_cons1) ."\t". sprintf("%.3f", $avg_rep1) ."\t".
		"$chr2\t$beg2\t$end2\t" . $ref->{strand_fosD} ."\t". $ref->{seq_fosD} ."\t". sprintf("%.3f", $avg_cons2) ."\t". sprintf("%.3f", $avg_rep2) ."\t".
		"$chr3\t$beg3\t$end3\t" . $ref->{strand_ETS} ."\t". $ref->{seq_ETS} ."\t". sprintf("%.3f", $avg_cons3) ."\t". sprintf("%.3f", $avg_rep3) ."\t".
		$ref->{begin_500} ."\t". $ref->{end_500} ."\t". $ref->{name} ."\t". 
		$ref->{gene_start_position} ."\t". $ref->{gene_end_position} ."\t". $ref->{gene_strand} ."\t".
		$ref->{location_type} ."\t". $ref->{gene_id} ."\t". $ref->{segment_id} . "\n";
}

exit(0);


### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "chip_seq_dmel";
        my $driver = "mysql";
        my $host = "andromeda.eb.local";
        my $username = "hox";
        my $password = "tfbs2006";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

### Connects to a database and returns databaseHandle
sub connect_to_db2
{
        my $databaseName = "mcginnes";
        my $driver = "mysql";
        my $host = "andromeda.eb.local";
        my $username = "hox";
        my $password = "tfbs2006";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh2 = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

