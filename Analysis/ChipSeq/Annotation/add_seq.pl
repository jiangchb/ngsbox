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
#  Module: Analysis::ChipSeq::Annotation::add_seq.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $dbh;
&connect_to_db();

my $dbh2;
&connect_to_db2();

### Load binding site table
my $q ="SELECT  chr,begin,end,strand,seq,grh_chr,begin_grh,end_grh,strand_grh,seq_grh,begin_500,end_500,segment_id,
		location_type,intron_position,segment_strand,gene_id,name,gene_start_position,gene_end_position,gene_strand
	#FROM   ETS2grh2segment2gene
	FROM	fosD2grh2segment2gene
	ORDER by chr, begin
	";
	my $sth = $dbh2->prepare($q);
	$sth->execute();
	
### Add average conservation and repeat probability for each TFBS
while(my $ref = $sth->fetchrow_hashref()) {
	my $chr1 = $ref->{'chr'};
	my $beg1 = $ref->{begin};
	my $end1 = $ref->{end};
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


        my $chr2 = $ref->{'grh_chr'};
        my $beg2 = $ref->{begin_grh};
        my $end2 = $ref->{end_grh};
        my $length2 = $end2 - $beg2 + 1;


        ### Load conservation, repeat status and sequence for second binding site
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

	print 	"$chr1\t$beg1\t$end1\t" . $ref->{strand} ."\t". $ref->{seq} ."\t". sprintf("%.3f", $avg_cons1) ."\t". sprintf("%.3f", $avg_rep1) ."\t". 
		"$chr2\t$beg2\t$end2\t" . $ref->{strand_grh} ."\t". $ref->{seq_grh} ."\t". sprintf("%.3f", $avg_cons2) ."\t". sprintf("%.3f", $avg_rep2) ."\t".
		$ref->{begin_500} ."\t". $ref->{end_500} ."\t". $ref->{segment_id} ."\t". $ref->{location_type} ."\t". 
		$ref->{intron_position} ."\t". $ref->{segment_strand} ."\t". $ref->{gene_id} ."\t". $ref->{name} ."\t". 
		$ref->{gene_start_position} ."\t". $ref->{gene_end_position} ."\t". $ref->{gene_strand} . "\n";
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

