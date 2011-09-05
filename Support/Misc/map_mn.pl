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
#  Module: Support::Misc::map_mn.pl
#  Purpose:
#  In:
#  Out:
#

####################################################################################
#Author 	Stephan Ossowski 
#Date 		05/15/07
#Version	0.9
#Input		Fragment file
#Function	Map a set of short sequence reads to a reference geneome allowing gaps
####################################################################################

use DBI;

my %seqs = ();
&read_fasta;

my $dbh;
&connect_to_db();

# Read each fragment from 2010 set
my $q = "SELECT distinct frag_id FROM mn_fragments
                WHERE ecotype = 'bur-0' && chromosome = 1
                ORDER BY id";
my $sth = $dbh->prepare($q);
$sth->execute();
while( my $ref = $sth->fetchrow_hashref() ) {
	my $frag_id        = $ref->{'frag_id'};
	my $frag_substring = $seqs{$frag_id};
	my $id             = $frag_id + 100000000;
	my $end_id         = $id + length($frag_substring);
	my $seq_ref        = "";

	# Read reference sequence
	my $q_col = "SELECT base FROM seq_ref WHERE id BETWEEN $id AND $end_id ORDER BY id";
	my $sth_col = $dbh->prepare($q_col);
	$sth_col->execute();
	while( my $ref_col = $sth_col->fetchrow_hashref() ) {
		$seq_ref .= $ref_col->{'base'};
	}

	# Write vmatch database file
	open MKV_FILE, ">$frag_id.db" or die "Could not open vmatch db file";
	print MKV_FILE ">$frag_id.db\n$seq_ref";
	close MKV_FILE;
	system("mkvtree -db $frag_id.db -dna -pl 1 -allout");


	# Write query file
	$frag_substring =~ s/^N*//;
	$frag_substring =~ s/N*$//;
	$frag_substring =~ s/-//g;
	open QUERY, ">$frag_id.fa" or die "Could not open query file";
	print QUERY ">$frag_id\n$frag_substring";
	close QUERY;
	
	#print STDERR "$frag_id\n$frag_substring\n";
		
	# Align reference and ecotype sequence using Vmatch
	my $length = length($frag_substring);
	my $error = int($length / 25);
	if($length > 60) { $error = 25; }
	if( ($frag_substring =~ /[ACTG]/) && ($length > 25) ) {
		system("vmatch -q $frag_id.fa -d -l $length -e $error -s -showdesc 30 -noevalue -noscore -noidentity $frag_id.db > $frag_id.align");
		&read_vmatch_output("$frag_id.align", $id);
	}
		
	# Clean up
	#system("rm $frag_id.align");
	system("rm $frag_id.fa");
	system("rm $frag_id.db*");
}
exit(0);


sub read_vmatch_output
{
	my $file = shift;
	my $id = shift;
	my $offset = 0;
	my $started = 0;
	open (MAP, "$file") or die;

	# Read command line
	my $junk = <MAP>;
	return(0) if not defined $junk;
	
	# Read statistics
	my $stats = <MAP>;
	return(0) if not defined $stats;

	chomp($stats);
	my @map_elem = split(/\s+/, $stats);
	my $total_length = $map_elem[0];
	my $chr = $map_elem[1];
	my $position = $map_elem[2] + 1;
	my $orientation = $map_elem[3];
	my $query_length = $map_elem[4];
	my $frag_name = $map_elem[5];
	my $vm_mismatches = $map_elem[7];
	$vm_mismatches = abs($vm_mismatches);

	# Read alignment
	while(<MAP>) {
		my $subject   = $_;
		if ( $subject eq "\n" ) { return(0); }
		my $polymorph = <MAP>;
		if( $polymorph =~ /Query/ ) { $junk = <MAP>; $offset += 60; next; }
		else {
			my $query = <MAP>;
			my ($junk1, $subject_seq, $junk2) = split(/\s/, $subject);
			my ($junk3, $query_seq, $junk4)   = split(/\s/, $query);
			my $pos = 0;
			for(my $x = 0; $x < length($subject_seq); $x++) {
				my $s = substr($subject_seq, $x, 1);
				my $q = substr($query_seq, $x, 1);
				if($s ne $q) { 
					if( ($q eq '-') && ($started == 1) ) {	# Deletion
						print "DELETION:  "; print $id + $pos + $offset; print "\t$frag_name\n";
					}
					elsif( ($s eq '-') && ($started == 1) ) {	# Insertion
						print "INSERTION: "; print $id + $pos + $offset; print "\t$frag_name\n";
						$pos--;
					}
					else {			# Gap
						if( ($s ne 'N') && ($q ne 'N') ) {
							print "SNP:       "; print $id + $pos + $offset; print "\t$frag_name\n";
						}
						$started = 1;
					}
				}
				$pos++;
			}
			$junk = <MAP>;
			$offset += 60;
		}
	}

	close(MAP);
	return(0);
}

sub read_fasta
{
        my $id = "";
        my @files = glob("*");
        foreach my $seqfile (@files) {
                open MN_IN, "$seqfile" or die "Could not open input file";
                while( <MN_IN> ) {
                        chomp;
                        if($_ =~ ">Bur-0") {
                                my $current_seq = <MN_IN>;
                                chomp($current_seq);
                                $seqs{$seqfile} .= $current_seq;
                        }
                }
                close(MN_IN);
        }

        return(0);
}

sub connect_to_db
{
        my $databaseName = "solexa";
        my $driver = "mysql";
        my $host = "localhost";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

