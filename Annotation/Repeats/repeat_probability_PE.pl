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
#  Module: Annotation::Repeats::repeat_probability_PE.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $database = shift;
my $file     = shift;

### Read insert size distribution ###
my $total_occurence = 0;
my $min_insert = 99999;
my $max_insert = 0;
my %probability = ();
open IN, $file or die "Cannot open input file\n";
while(<IN>) {
	chomp;
	my $line = $_;
	$line =~ s/^\s+//g;
	my ($occurence, $pos) = split(/ /, $line);
	$probability{$pos} = $occurence;
	$total_occurence += $occurence;
	if($pos < $min_insert) { $min_insert = $pos; }
	if($pos > $max_insert) { $max_insert = $pos; }
}
foreach(sort {$a <=> $b} keys %probability) {
	$probability{$_} = $probability{$_} / $total_occurence;
}
close(IN);

my $dbh='';
&connect_to_db();
my $q = "SELECT chromosome, max_pos FROM seq_max";
my $sth = $dbh->prepare($q);
$sth->execute();

### For each chromosome of reference genome
while( my $ref = $sth->fetchrow_hashref() ) {
	my $chr     = $ref->{chromosome};
	my $max_pos = $ref->{max_pos};

	### For all position in single end repeat table
	for (my $win = 251; $win <= $max_pos; $win += 10000) {
		my $current_end = $win + 10000 - 1;
		if($current_end > $max_pos) { $current_end = $max_pos; }

		my %cov_hash = ();
		my %hit_hash = ();
	
		my $q2="SELECT 	position, average_hits, coverage
			FROM 	repeat_seg_36
			WHERE 	chromosome = $chr &&
				position BETWEEN $win -250 and $current_end + 250
			ORDER BY position";
		my $sth2 = $dbh->prepare($q2);
		$sth2->execute();
	
		while( my $ref2 = $sth2->fetchrow_hashref() ) {
			$cov_hash{$ref2->{position}} = $ref2->{coverage};
			$hit_hash{$ref2->{position}} = $ref2->{average_hits};
		}

		### Calculate repeat rescue probability
		for(my $i = $win; $i <= $current_end; $i++) {
			my $sum_repeat_up = 0;
			my $sum_repeat_down = 0;
			my $pos_count = $max_insert;
			if(defined $hit_hash{$i} && $hit_hash{$i} > 1) {

				for(my $n = $i - $max_insert; $n <= $i - $min_insert; $n++) {
					if(defined $hit_hash{$n} && $hit_hash{$n} == 1) {
						$sum_repeat_up += $probability{$pos_count};
					}
					$pos_count--;
				}

				$pos_count = $min_insert;
				for(my $n = $i + $min_insert; $n <= $i+$max_insert; $n++) {
					if(defined $hit_hash{$n} && $hit_hash{$n} == 1) {
						$sum_repeat_down += $probability{$pos_count};
					}
					$pos_count++;
				}

				my $result = ($sum_repeat_up + $sum_repeat_down) / 2;
				print "$chr\t$i\t$result\n";
			}
		}
	}
}

#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
	my $databaseName = $database;
        my $driver = "mysql";
	my $host = "ume.fml.local";
	my $username = "solexa";
	my $password = "s0lexa";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}


