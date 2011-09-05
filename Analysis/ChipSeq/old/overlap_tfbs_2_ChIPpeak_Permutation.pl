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
#  Module: Analysis::ChipSeq::old::overlap_tfbs_2_ChIPpeak_Permutation.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $peak           = shift;
my $enriched_table = shift;
my $tfbs_in        = shift;

my @tfbs_tables = split(/,/, $tfbs_in);

my %total_good = ();
my %total_bad  = ();
my $dbh;
&connect_to_db();

### Query chromosome length
my $q0 = "SELECT chromosome, chr_alias, max_pos FROM seq_max where chromosome BETWEEN 3 AND 8 ORDER BY chromosome";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();

### For each chromosome
while( my $ref0 = $sth0->fetchrow_hashref() ) {
	my $chr       = $ref0->{chromosome};
	my $chr_alias = $ref0->{chr_alias};
	my $max_pos   = $ref0->{max_pos};
	my %tfbs = ();

	### Load TFBS prediction
	foreach my $tfbs_table (@tfbs_tables) {
                my $q= "SELECT  begin, end
                        FROM    $tfbs_table
                        WHERE   chromosome = '$chr_alias'
                        ORDER BY begin
        	";
        	my $sth = $dbh->prepare($q);
        	$sth->execute();

        	while(my $ref = $sth->fetchrow_hashref()) {
        	        for(my $i = $ref->{begin} - 100; $i <= $ref->{end} + 100; $i++) {
        	                $tfbs{$i} = 1;
        	        }
        	}
		$sth->finish();
	}

        ### Load ChipSeq peaks
        my $q=" SELECT distinct begin, end
                FROM    $enriched_table
                WHERE   chromosome = $chr &&
                        max_cov >= $peak &&
                        segment_length >= 100 &&
			fbgn in (select sequence_derived_from from array_dfd_anno_gene_jun08 where dfd_early_no_rec_hs > 1.5 || dfd_early_no_rec_hs between 0.001 and 0.5 || dfd_early_gof > 1.5 || dfd_early_gof between 0.001 and 0.5 || dfd_late_gof > 1.5 || dfd_late_gof between 0.001 and 0.5)
                ORDER by chromosome, begin
        ";
        my $sth = $dbh->prepare($q);
        $sth->execute();

	my %peak = ();
	while(my $ref = $sth->fetchrow_hashref()) {
		my $beg = $ref->{begin};
		my $end = $ref->{end};
		$peak{$beg} = $end;
	}
	$sth->finish();

	### Compare
	for( my $permute = 1; $permute < 1001; $permute++) {
		my %result = (good => 0, bad => 0);
		my $curr_rand = int(rand($max_pos)) + 1;
		
		foreach my $beg (keys %peak) {
			my $end = $peak{$beg};
			my $rand_beg = (($beg + $curr_rand)%$max_pos) + 1;
			my $rand_end = (($end + $curr_rand)%$max_pos) + 1;
			my $hit = 0;

			for(my $pos = $rand_beg; $pos <= $rand_end; $pos++) {
				if(exists $peak{$pos}) { $hit = 1; }
			}
			if($hit == 1) { $result{good}++; }
			else          { $result{bad}++; }
		}
		$total_good{$permute} += $result{good};
		$total_bad{$permute} += $result{bad};
	}
}
$sth0->finish();
	
foreach (sort keys %total_good) {
	print "$total_good{$_}\t$total_bad{$_}\t". sprintf("%.2f", $total_good{$_}/$total_bad{$_}*100) ."%\n";
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

