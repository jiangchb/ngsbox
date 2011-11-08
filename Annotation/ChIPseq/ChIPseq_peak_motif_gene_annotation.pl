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
#  Module: Annotation::ChIPseq::ChIPseq_peak_motif_gene_annotation.pl
#  Purpose:
#  In:
#  Out:
#

use POSIX;
use DBI;

my $usage = "\n\n$0 database\n\n";
my $database       = shift or die $usage;

my $dbh = '';
&connect_to_db();


### Read gene annotation
my $ann_enst2geneID = $dbh->selectall_hashref('SELECT * FROM ann_enst2geneID', 'enst_name');
#foreach my $id (keys %$ann_enst2geneID) {
#	print "Value of ID $id is $ann_enst2geneID->{$id}->{ensg_name}\n";
#}



### Read TSS annotation
my $ann_switchGearTSS_fwd = $dbh->selectall_hashref('SELECT * FROM ann_switchGearTSS WHERE strand = "+"', 'name');
my $ann_switchGearTSS_rev = $dbh->selectall_hashref('SELECT * FROM ann_switchGearTSS WHERE strand = "-"', 'name');
#foreach my $id (keys %$ann_switchGearTSS_fwd) {
#	print "Value of ID $id is $ann_switchGearTSS_fwd->{$id}->{strand}\n";
#}



### Read fimo motif prediction
my $q= "SELECT * from fimo_gw";
my $sth = $dbh->prepare($q);
$sth->execute();
my %fimo = ();
while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr_fimo};
	my $beg = $ref->{start_fimo};
	my $end = $ref->{stop_fimo};
	$fimo{"$chr#$beg#$end"} = \%$ref;
	#print $fimo{"$chr#$beg#$end"}->{seq_fimo} . "\n";
}

#foreach my $id (keys %fimo) {
#	print "$id\t" . $fimo{$id}->{seq_fimo} . "\n";
#}



### Read BED file
$q= "SELECT * from peak";
$sth = $dbh->prepare($q);
$sth->execute();
my %peak = ();
while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr};
	my $beg = $ref->{start};
	my $end = $ref->{end};
	$peak{"$chr#$beg#$end"} = \%$ref;
}

### Read PeakAnalyzer NDG
$q= "SELECT * from pa_ndg";
$sth = $dbh->prepare($q);
$sth->execute();
my %pa_ndg = ();
while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr_pa};
	my $beg = $ref->{start_pa};
	my $end = $ref->{end_pa};
	my $peak_middle = floor( ($beg + $end) / 2 );
	
	my $enst_down_fwd = $ref->{Downstream_FW_Gene};
	my $enst_down_rev = $ref->{Downstream_REV_Gene};


	### Get fimo motif
	my $fimo_q= "SELECT * FROM fimo_gw WHERE score_fimo >= 6 && chr_fimo = '$chr' && start_fimo between $beg and $end";
	my $fimo_sth = $dbh->prepare($fimo_q);
	$fimo_sth->execute();
	my %fimo = ();
	my $TFBS_string = "";
	while (my $ref = $fimo_sth->fetchrow_hashref()) {
		my $peak_middle_dist = $peak_middle - $ref->{start_fimo};

		$TFBS_string .= "(". $ref->{seq_fimo} .",$peak_middle_dist,". $ref->{score_fimo}. "),";
	}
	chop($TFBS_string);
	if($TFBS_string eq "") { $TFBS_string = "NA"; }
	

	### TODO if motif found, than check TSS relative to motif middle, otherwise relative to peak middle


	### Get downstream fwd TSS
	my $TSS_name_fwd = "";
	my $TSS_dist_fwd = -1;
	my $TSS_name_rev = "";
	my $TSS_dist_rev = -1;

	my $ann_switchGearTSS_fwd = $dbh->selectall_hashref("SELECT * FROM ann_switchGearTSS WHERE strand = '+' AND chrom = '$chr' AND chromStart >= $beg ORDER by chromStart ASC limit 1", "name");
	foreach my $id (keys %$ann_switchGearTSS_fwd) {
		$TSS_name_fwd = $ann_switchGearTSS_fwd->{$id}->{name};
		$TSS_dist_fwd = $ann_switchGearTSS_fwd->{$id}->{chromStart} - $beg; 
	}

	my $ann_switchGearTSS_rev = $dbh->selectall_hashref("SELECT * FROM ann_switchGearTSS WHERE strand = '-' AND chrom = '$chr' AND chromStart <= $end ORDER by chromStart DESC limit 1", "name");
	foreach my $id (keys %$ann_switchGearTSS_rev) { 
		$TSS_name_rev = $ann_switchGearTSS_rev->{$id}->{name};
		$TSS_dist_rev = $end - $ann_switchGearTSS_rev->{$id}->{chromStart}; 
	}



	### Print peak features
	print   $ref->{venn_set} ."\t". $ref->{venn_subset} ."\t".
		$ref->{chr_pa} ."\t". $ref->{start_pa} ."\t". $ref->{end_pa} ."\t".
		$peak{"$chr#$beg#$end"}->{region} ."\t". $peak{"$chr#$beg#$end"}->{peakscore} ."\t".
		$ref->{overlaped_transcripts} ."\t";

	### Print fwd TSS annotation
	print	"$TSS_name_fwd\t$TSS_dist_fwd\t";

	### Print fwd gene annotation
	if(exists $ann_enst2geneID->{$enst_down_fwd}->{ensg_name}) {
		print	$enst_down_fwd ."\t". $ann_enst2geneID->{$enst_down_fwd}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{GeneID} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{refseq} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{symbol} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{name_desc} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{description} ."\t";
	}
	else {
		print	$enst_down_fwd ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";
	}

	### Print rev TSS annotation
	print	"$TSS_name_rev\t$TSS_dist_rev\t";

	### Print rev gene annotation
	if(exists $ann_enst2geneID->{$enst_down_rev}->{ensg_name}) {
		print	$enst_down_rev ."\t". $ann_enst2geneID->{$enst_down_rev}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{GeneID} ."\t". $ann_enst2geneID->{$enst_down_rev}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{refseq} ."\t". $ann_enst2geneID->{$enst_down_rev}->{symbol} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{name_desc} ."\t". $ann_enst2geneID->{$enst_down_rev}->{description} ."\t";
	}
	else {
		print	$enst_down_rev ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";
	}

	### Print motif information
	print "$TFBS_string\n";

}

exit(0);



### Read PeakAnalyzer Overlap
$q= "SELECT * from pa_overlap";
$sth = $dbh->prepare($q);
$sth->execute();
my %pa_ovl = ();
while (my $ref = $sth->fetchrow_hashref()) {

	my $chr = $ref->{chr_pa};
	my $beg = $ref->{start_pa};
	my $end = $ref->{end_pa};
	# $pa_ovl{"$chr#$beg#$end"} = \%$ref;
	
	my $enst = $ref->{overlaped_transcripts};

	if(exists $ann_enst2geneID->{$enst}->{ensg_name}) {
	
		print 	$ref->{venn_set} ."\t". $ref->{venn_subset} ."\t". 
			$ref->{chr_pa} ."\t". $ref->{start_pa} ."\t". $ref->{end_pa} ."\t".
			$peak{"$chr#$beg#$end"}->{region} ."\t". $peak{"$chr#$beg#$end"}->{peakscore} ."\t".
			$enst ."\t". $ann_enst2geneID->{$enst}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst}->{GeneID} ."\t". $ann_enst2geneID->{$enst}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst}->{refseq} ."\t". $ann_enst2geneID->{$enst}->{symbol} ."\t".
			$ann_enst2geneID->{$enst}->{name_desc} ."\t". $ann_enst2geneID->{$enst}->{description} ."\t".
			$ref->{Overlap_Begin} ."\t". $ref->{Overlap_Center} ."\t".$ref->{Overlap_End} ."\n";
	}
	else {
		print   $ref->{venn_set} ."\t". $ref->{venn_subset} ."\t".
			$ref->{chr_pa} ."\t". $ref->{start_pa} ."\t". $ref->{end_pa} ."\t".
			$ref->{overlaped_transcripts} ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t".
			$ref->{Overlap_Begin} ."\t". $ref->{Overlap_Center} ."\t".$ref->{Overlap_End} ."\n";
	}
}
#foreach my $id (keys %pa_ovl) {
#	print "$id\t" . $pa_ovl{$id}->{overlaped_genes} . "\n";
#}

exit(0);



#####################################################
### Connects to a database and returns databaseHandle
#####################################################
sub connect_to_db
{
	my $databaseName = "$database";
        my $driver = "mysql";
	my $host = "pou.crg.es";
	my $username = "gdvisitor";
	my $password = "";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to genopheno-db. Connect error: $DBI::errstr";
}


