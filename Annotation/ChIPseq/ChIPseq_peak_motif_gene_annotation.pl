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



### Read BED file
my $q= "SELECT * from peak";
my $sth = $dbh->prepare($q);
$sth->execute();
my %peak = ();
while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr};
	my $beg = $ref->{start} - 1;
	my $end = $ref->{end} - 1;
	$peak{"$chr#$beg#$end"} = \%$ref;
}



### Peak annotation (using PeakAnalyzer NDG, Homer and UCSC hg19)
$q= "SELECT * from pa_ndg";
$sth = $dbh->prepare($q);
$sth->execute();
my %pa_ndg = ();

open NDG, ">peak_annotation_pa_ndg.txt" or die "Cannot open output file";
print NDG "Set\tSubset\tChr\tPeak_Start\tPeak_End\tPeak_Region\tPeak_Score\tOverlapped_Transcripts\t".
        "Homer_Annotation\tHomer_Nearest_PromoterID\tHomer_Nearest_Promoter_Distance\tHomer_Nearest_Promoter_EntrezID\ttHomer_Nearest_Promoter_ENSG\tHomer_Nearest_Promoter_Symbol\t".
        "SwitchGear_TSS_FWD\tSwitchGear_TSS_Distance_FWD\tSwitchGear_TSS_REV\tSwitchGear_TSS_Distance_REV\t".
        "ENST_FWD\tENSG_FWD\tEntrezID_FWD\tUniGene_FWD\tRefSeq_FWD\tSymbol_FWD\tName_FWD\tDescription_FWD\t".
        "ENST_REV\tENSG_REV\tEntrezID_REV\tUniGene_REV\tRefSeq_REV\tSymbol_REV\tName_REV\tDescription_REV\t".
        "Motifs(sequence,distance_from_peak,fimo_score)\n";

while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr_pa};
	my $beg = $ref->{start_pa} - 1;
	my $end = $ref->{end_pa} - 1;
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
	


	### Homer style annotation of feature type
	my $homer_q = "SELECT * FROM ann_homer_basic WHERE chr = 'chr$chr' AND ($beg BETWEEN start AND end || $end BETWEEN start AND end || start BETWEEN $beg AND $end)";
	my $homer_sth = $dbh->prepare($homer_q);
	$homer_sth->execute();
	my %homer = ();
	my $homer_string = "";
	while (my $ref = $homer_sth->fetchrow_hashref()) {
		$homer{$ref->{type}} .= $ref->{anno} . ",";
	}
	if(exists $homer{"P"})			{$homer_string = $homer{"P"}}
	elsif(exists $homer{"E"})		{$homer_string = $homer{"E"}}
	elsif(exists $homer{"5UTR"})		{$homer_string = $homer{"5UTR"}}
	elsif(exists $homer{"3UTR"})		{$homer_string = $homer{"3UTR"}}
	elsif(exists $homer{"I"})		{$homer_string = $homer{"I"}}
	elsif(exists $homer{"CpG-Island"})	{$homer_string = $homer{"CpG-Island"}}
	else{ $homer_string = "Intergenic,"; }
	chop($homer_string);


	### Get nearest PromoterID from homer hg19 (from all up/down, fwd/rev)
	my $promoter_name = "NA";
	my $promoter_id = "NA";
	my $promoter_symbol = "NA";
	my $promoter_ensg = "NA";
	my $promoter_dist = 999999999;


	my $promoter_up = $dbh->selectall_hashref("SELECT * FROM ann_homer_basic WHERE type = 'P' AND chr = 'chr$chr' AND start >= $peak_middle ORDER by start ASC limit 1", "anno");
	foreach my $id (keys %$promoter_up) {
		my $promoter_middle = floor( ($promoter_up->{$id}->{start} + $promoter_up->{$id}->{end}) / 2 );
		my $current_dist = abs($promoter_middle - $peak_middle);
		if($current_dist < $promoter_dist) {
			$promoter_dist = $current_dist;
			$promoter_name = $promoter_up->{$id}->{anno};
		}
	}
	my $promoter_down = $dbh->selectall_hashref("SELECT * FROM ann_homer_basic WHERE type = 'P' AND chr = 'chr$chr' AND start <= $peak_middle ORDER by start DESC limit 1", "anno");
	foreach my $id (keys %$promoter_down) {
		my $promoter_middle = floor( ($promoter_down->{$id}->{start} + $promoter_down->{$id}->{end}) / 2 );
		my $current_dist = abs($promoter_middle - $peak_middle);
		if($current_dist < $promoter_dist) {
			$promoter_dist = $current_dist;
			$promoter_name = $promoter_down->{$id}->{anno};
		}
	}
	$promoter_name =~ s/promoter-TSS \(//g;
	$promoter_name =~ s/\)//g;

	my $translate = $dbh->selectall_hashref("SELECT * FROM ann_enst2geneID where refseq = '$promoter_name' limit 1", "enst_name");
	foreach my $id (keys %$translate) {
		$promoter_id = $translate->{$id}->{GeneID};
		$promoter_symbol = $translate->{$id}->{symbol};
		$promoter_ensg = $translate->{$id}->{ensg_name};
	}
	

	### Get nearest TSS from SwitchGear (fwd and rev)
	my $TSS_name_fwd = "";
	my $TSS_dist_fwd = -1;
	my $TSS_name_rev = "";
	my $TSS_dist_rev = -1;

	my $ann_switchGearTSS_fwd = $dbh->selectall_hashref("SELECT * FROM ann_switchGearTSS WHERE strand = '+' AND chrom = '$chr' AND chromStart >= $beg ORDER by chromStart ASC limit 1", "name");
	foreach my $id (keys %$ann_switchGearTSS_fwd) {
		$TSS_name_fwd = $ann_switchGearTSS_fwd->{$id}->{name};
		$TSS_dist_fwd = $ann_switchGearTSS_fwd->{$id}->{chromStart} - $peak_middle; 
	}

	my $ann_switchGearTSS_rev = $dbh->selectall_hashref("SELECT * FROM ann_switchGearTSS WHERE strand = '-' AND chrom = '$chr' AND chromStart <= $end ORDER by chromStart DESC limit 1", "name");
	foreach my $id (keys %$ann_switchGearTSS_rev) { 
		$TSS_name_rev = $ann_switchGearTSS_rev->{$id}->{name};
		$TSS_dist_rev = $peak_middle - $ann_switchGearTSS_rev->{$id}->{chromStart}; 
	}


	### Print peak features
	print NDG $ref->{venn_set} ."\t". $ref->{venn_subset} ."\t".
		"chr$chr\t$beg\t$end\t".
		$peak{"$chr#$beg#$end"}->{region} ."\t". $peak{"$chr#$beg#$end"}->{peakscore} ."\t".
		$ref->{overlaped_transcripts} . "\t".
		"$homer_string\t$promoter_name\t$promoter_dist\t$promoter_id\t$promoter_ensg\t$promoter_symbol\t";

	### Print SwitchGear TSS annotation
	print NDG "$TSS_name_fwd\t$TSS_dist_fwd\t$TSS_name_rev\t$TSS_dist_rev\t";

	### Print fwd gene annotation
	if(exists $ann_enst2geneID->{$enst_down_fwd}->{ensg_name}) {
		print NDG $enst_down_fwd ."\t". $ann_enst2geneID->{$enst_down_fwd}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{GeneID} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{refseq} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{symbol} ."\t".
			$ann_enst2geneID->{$enst_down_fwd}->{name_desc} ."\t". $ann_enst2geneID->{$enst_down_fwd}->{description} ."\t";
	}
	else {
		print NDG $enst_down_fwd ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";
	}

	### Print rev gene annotation
	if(exists $ann_enst2geneID->{$enst_down_rev}->{ensg_name}) {
		print NDG $enst_down_rev ."\t". $ann_enst2geneID->{$enst_down_rev}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{GeneID} ."\t". $ann_enst2geneID->{$enst_down_rev}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{refseq} ."\t". $ann_enst2geneID->{$enst_down_rev}->{symbol} ."\t".
			$ann_enst2geneID->{$enst_down_rev}->{name_desc} ."\t". $ann_enst2geneID->{$enst_down_rev}->{description} ."\t";
	}
	else {
		print NDG $enst_down_rev ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";
	}

	### Print motif information
	print NDG "$TFBS_string\n";

}



### Read PeakAnalyzer Overlap
$q= "SELECT * from pa_overlap";
$sth = $dbh->prepare($q);
$sth->execute();
my %pa_ovl = ();
open OVL, ">peak_annotation_pa_ovl.txt" or die "Cannot open output file";
while (my $ref = $sth->fetchrow_hashref()) {

	my $chr = $ref->{chr_pa};
	my $beg = $ref->{start_pa} - 1;
	my $end = $ref->{end_pa} - 1;
	# $pa_ovl{"$chr#$beg#$end"} = \%$ref;
	
	my $enst = $ref->{overlaped_transcripts};

	if(exists $ann_enst2geneID->{$enst}->{ensg_name}) {
	
		print OVL $ref->{venn_set} ."\t". $ref->{venn_subset} ."\t". 
			$ref->{chr_pa} ."\t". $ref->{start_pa} ."\t". $ref->{end_pa} ."\t".
			$peak{"$chr#$beg#$end"}->{region} ."\t". $peak{"$chr#$beg#$end"}->{peakscore} ."\t".
			$enst ."\t". $ann_enst2geneID->{$enst}->{ensg_name} ."\t".
			$ann_enst2geneID->{$enst}->{GeneID} ."\t". $ann_enst2geneID->{$enst}->{Unigene} ."\t".
			$ann_enst2geneID->{$enst}->{refseq} ."\t". $ann_enst2geneID->{$enst}->{symbol} ."\t".
			$ann_enst2geneID->{$enst}->{name_desc} ."\t". $ann_enst2geneID->{$enst}->{description} ."\t".
			$ref->{Overlap_Begin} ."\t". $ref->{Overlap_Center} ."\t".$ref->{Overlap_End} ."\n";
	}
	else {
		print OVL $ref->{venn_set} ."\t". $ref->{venn_subset} ."\t".
			$ref->{chr_pa} ."\t". $ref->{start_pa} ."\t". $ref->{end_pa} ."\t".
			$ref->{overlaped_transcripts} ."\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t".
			$ref->{Overlap_Begin} ."\t". $ref->{Overlap_Center} ."\t".$ref->{Overlap_End} ."\n";
	}
}

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

