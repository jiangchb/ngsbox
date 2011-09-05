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
#  Module: Analysis::ChipSeq::Annotation::average_conservation.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $dbh;
&connect_to_db();

### Load average conservation
my $q ="SELECT  id, parent_cluster_id, type, begin, end, chromosome, strand, average_conservation, patser_score, prediction_tool
	FROM    tfbs_cisA_hoxB
	ORDER by chromosome, begin
	";
	my $sth = $dbh->prepare($q);
	$sth->execute();
	
### Add average conservation and repeat probability for each TFBS
while(my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chromosome};
	my $beg = $ref->{begin};
	my $end = $ref->{end};
	my $length = $end - $beg + 1;

	my $patser_score = "\\N";
	if( defined $ref->{patser_score} ) {
		$patser_score = $ref->{patser_score};
	}

	### Load conservation, repeat status and sequence from seq_ref
        my $q2="SELECT  conservation, isRepeat, nucleotide
                FROM    seq_ref
                WHERE   chromosome = '$chr' AND
			position BETWEEN $beg AND $end
                ORDER BY chromosome, position
        ";
        my $sth2 = $dbh->prepare($q2);
        $sth2->execute();

	my $cons_sum = 0;
	my $rep_sum = 0;
	my $seq = "";
	while(my $ref2 = $sth2->fetchrow_hashref()) {
		if(defined $ref2->{conservation}) {
			$cons_sum += $ref2->{conservation};
		}
		else { 
			$cons_sum += 0;
		}
		$rep_sum += $ref2->{isRepeat};
		$seq .= $ref2->{nucleotide};
	}
	my $avg_cons = $cons_sum / $length;
	my $avg_rep  = $rep_sum / $length;

	print 	$ref->{id} ."\t". $ref->{parent_cluster_id} ."\t". $ref->{type} ."\t$beg\t$end\t$chr\t" . 
		$ref->{strand} . "\t$avg_cons\t$patser_score\t". $ref->{prediction_tool} . "\t$avg_rep\t$seq\n";
}

exit(0);


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

