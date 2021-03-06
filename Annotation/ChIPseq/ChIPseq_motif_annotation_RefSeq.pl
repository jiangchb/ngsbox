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
my $database = shift or die $usage;

my $dbh = '';
&connect_to_db();


### Read gene annotation translation
my $ann_enst2geneID       = $dbh->selectall_hashref('SELECT * FROM ann_enst2geneID group by refseq', 'refseq');
my $ann_NCBI_gene2ensembl = $dbh->selectall_hashref('SELECT * FROM ann_NCBI_gene2ensembl', 'RNA_nucleotide_accession');
my $ann_NCBI_gene_info    = $dbh->selectall_hashref('SELECT * FROM ann_NCBI_gene_info', 'GeneID');




### Print peak annotation header
open NDG, ">motif_annotation_refseq.txt" or die "Cannot open output file";
print NDG "Chr\tMotif_Start\tMotif_End\tMotif_Sequence\tMotif_Score\t".

	# Homer style
	"Homer:Annotation\tHomer:Nearest_Promoter_Distance\tHomer:Nearest_Promoter_RefSeq\t".
	"Homer:Nearest_Promoter_GeneID\tHomer:Nearest_Promoter_EnsG\t".

	# SwitchGearTSS
	"SwitchGear:TSS_FWD\tSwitchGear:TSS_Distance_FWD\tSwitchGear:TSS_REV\tSwitchGear:TSS_Distance_REV\t".

	# RefSeqGenes
	"RefSeq:Overlap\t".
	"RefSeq:FWD\tDistance:txStart\tDistance:cdsStart\tSymbol:FWD\tGeneID:FWD\tEnsT:FWD\tEnsG:FWD\tSyonyms:FWD\tDescription:FWD\t".
	"RefSeq:REV\tDistance:txStart\tDistance:cdsStart\tSymbol:REV\tGeneID:REV\tEnsT:REV\tEnsG:REV\tSyonyms:REV\tDescription:REV\n";


### Get motifs
my $q= "SELECT * from fimo_gw WHERE score_fimo >= 6";
my $sth = $dbh->prepare($q);
$sth->execute();
my %peak = ();

while (my $ref = $sth->fetchrow_hashref()) {
	my $chr = $ref->{chr_fimo};
	my $beg = $ref->{start_fimo};
	my $end = $ref->{stop_fimo};
	my $peak_middle = ceil( ($beg + $end) / 2 );



	### Get RefSeqGenes fwd downstream gene
	my $refseq_fwd = "NA";
	my $ann_UCSC_RefSeqGenes_fwd = $dbh->selectall_hashref("SELECT * FROM ann_UCSC_RefSeqGenes WHERE strand = '+' AND chrom = 'chr$chr' AND txStart >= $beg ORDER by txStart ASC limit 1", "name");
	foreach my $id (keys %$ann_UCSC_RefSeqGenes_fwd) {
		$refseq_fwd = $ann_UCSC_RefSeqGenes_fwd->{$id}->{name};
	}
	my $refseq_fwd_dist_txStart  = $ann_UCSC_RefSeqGenes_fwd->{$refseq_fwd}->{txStart} - $peak_middle;
	my $refseq_fwd_dist_cdsStart = $ann_UCSC_RefSeqGenes_fwd->{$refseq_fwd}->{cdsStart} - $peak_middle;


	### Get RefSeqGenes rev downstream gene
	my $refseq_rev = "NA";
	my $ann_UCSC_RefSeqGenes_rev = $dbh->selectall_hashref("SELECT * FROM ann_UCSC_RefSeqGenes WHERE strand = '-' AND chrom = 'chr$chr' AND txEnd <= $end ORDER by txEnd DESC limit 1", "name");
	foreach my $id (keys %$ann_UCSC_RefSeqGenes_rev) {
		$refseq_rev = $ann_UCSC_RefSeqGenes_rev->{$id}->{name};
	}
	my $refseq_rev_dist_txStart  = $peak_middle - $ann_UCSC_RefSeqGenes_rev->{$refseq_rev}->{txEnd};
	my $refseq_rev_dist_cdsStart = $peak_middle - $ann_UCSC_RefSeqGenes_rev->{$refseq_rev}->{cdsEnd};

	
	### Get RefSeqGenes overlapping genes
	my $refseq_ovl = "";
	my %refseq_ovl_hash = ();
	my $ann_UCSC_RefSeqGenes_ovl = $dbh->selectall_hashref("SELECT * FROM ann_UCSC_RefSeqGenes WHERE chrom = 'chr$chr' AND $peak_middle between cdsStart and cdsEnd", "name");
	foreach my $id (keys %$ann_UCSC_RefSeqGenes_ovl) {
		$refseq_ovl_hash{$ann_UCSC_RefSeqGenes_ovl->{$id}->{name}} = 1;
	}
	foreach my $refseq (sort keys %refseq_ovl_hash) {
		$refseq_ovl .= "$refseq,";
	}
	chop($refseq_ovl);
	if($refseq_ovl eq "") { $refseq_ovl = "NA"; }



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
	my $promoter_ensg = "NA";
	my $promoter_dist = 999999999;


	my $promoter_up = $dbh->selectall_hashref("SELECT * FROM ann_homer_hg19_tss WHERE chr = 'chr$chr' AND start + 2000 >= $peak_middle ORDER by start ASC limit 1", "refseq");
	foreach my $id (keys %$promoter_up) {
		my $promoter_middle = $promoter_up->{$id}->{start} + 2001;
		my $current_dist = $promoter_middle - $peak_middle;
		if($promoter_up->{$id}->{strand} eq "0") {
			$current_dist = $peak_middle - $promoter_middle;
		}
		if(abs($current_dist) < abs($promoter_dist)) {
			$promoter_dist = $current_dist;
			$promoter_name = $promoter_up->{$id}->{refseq};
		}
	}
	my $promoter_down = $dbh->selectall_hashref("SELECT * FROM ann_homer_hg19_tss WHERE chr = 'chr$chr' AND start + 2000 <= $peak_middle ORDER by start DESC limit 1", "refseq");
	foreach my $id (keys %$promoter_down) {
		my $promoter_middle = $promoter_down->{$id}->{start} + 2001;
		my $current_dist = $promoter_middle - $peak_middle;
		if($promoter_down->{$id}->{strand} eq "0") {
			$current_dist = $peak_middle - $promoter_middle;
		}
		if(abs($current_dist) < abs($promoter_dist)) {
			$promoter_dist = $current_dist;
			$promoter_name = $promoter_down->{$id}->{refseq};
		}
	}

	my $translate = $dbh->selectall_hashref("SELECT * FROM ann_homer_gene where refseq = '$promoter_name' || symbol = '$promoter_name' limit 1", "symbol");
	foreach my $id (keys %$translate) {
		$promoter_id = $translate->{$id}->{GeneID};
		$promoter_ensg = $translate->{$id}->{ensg_name};
		$promoter_name = $translate->{$id}->{refseq};
	}

	

	### Get nearest TSS from SwitchGear (fwd and rev)
	my $TSS_name_fwd = "";
	my $TSS_dist_fwd = 999999999;
	my $TSS_name_rev = "";
	my $TSS_dist_rev = 999999999;

	my $ann_UCSC_switchGearTSS_fwd = $dbh->selectall_hashref("SELECT * FROM ann_UCSC_switchGearTSS WHERE isPseudo = 0 AND confScore > 5 AND strand = '+' AND chrom = '$chr' AND chromStart >= $beg ORDER by chromStart ASC limit 5", "name");
	foreach my $id (keys %$ann_UCSC_switchGearTSS_fwd) {
		my $current_dist = $ann_UCSC_switchGearTSS_fwd->{$id}->{chromStart} - $peak_middle;
		if( abs($current_dist) < abs($TSS_dist_fwd) ) {
			$TSS_name_fwd = $ann_UCSC_switchGearTSS_fwd->{$id}->{name};
			$TSS_dist_fwd = $current_dist;
		}
	}

	my $ann_UCSC_switchGearTSS_rev = $dbh->selectall_hashref("SELECT * FROM ann_UCSC_switchGearTSS WHERE isPseudo = 0 AND confScore > 5 AND strand = '-' AND chrom = '$chr' AND chromStart <= $end ORDER by chromStart DESC limit 5", "name");
	foreach my $id (keys %$ann_UCSC_switchGearTSS_rev) {
		my $current_dist = $peak_middle - $ann_UCSC_switchGearTSS_rev->{$id}->{chromStart};
		if( abs($current_dist) < abs($TSS_dist_rev) ) {
			$TSS_name_rev = $ann_UCSC_switchGearTSS_rev->{$id}->{name};
			$TSS_dist_rev = $current_dist;
		}
	}




	### Print peak features
	print NDG "chr$chr\t$beg\t$end\t" . $ref->{seq_fimo} ."\t". $ref->{score_fimo} ."\t";

	### Print Homer style annotation
	print NDG "$homer_string\t$promoter_dist\t$promoter_name\t$promoter_id\t$promoter_ensg\t";

	### Print SwitchGear TSS annotation
	print NDG "$TSS_name_fwd\t$TSS_dist_fwd\t$TSS_name_rev\t$TSS_dist_rev\t";

	### Print overlapping genes (CDS & Intron)
	print NDG "$refseq_ovl\t";

	### Print fwd gene annotation
	print NDG "$refseq_fwd\t$refseq_fwd_dist_txStart\t$refseq_fwd_dist_cdsStart\t". 
			$ann_UCSC_RefSeqGenes_fwd->{$refseq_fwd}->{name2} ."\t";

	if( ($refseq_fwd ne "NA") && (exists $ann_NCBI_gene2ensembl->{$refseq_fwd}->{GeneID}) ) {
		my $geneID = $ann_NCBI_gene2ensembl->{$refseq_fwd}->{GeneID};
		
		print NDG $geneID ."\t".
			$ann_NCBI_gene2ensembl->{$refseq_fwd}->{Ensembl_rna_id} ."\t".
			$ann_NCBI_gene2ensembl->{$refseq_fwd}->{Ensembl_gene_id} ."\t".
			$ann_NCBI_gene_info->{$geneID}->{Synonyms} ."\t".
			$ann_NCBI_gene_info->{$geneID}->{description} ."\t";
	}

	elsif( ($refseq_fwd ne "NA") && (exists $ann_enst2geneID->{$refseq_fwd}->{GeneID}) ) {
		my $geneID = $ann_enst2geneID->{$refseq_fwd}->{GeneID};
		
		print NDG $geneID ."\t".
			$ann_enst2geneID->{$refseq_fwd}->{enst_name} ."\t".
			$ann_enst2geneID->{$refseq_fwd}->{ensg_name} ."\t".
			$ann_enst2geneID->{$refseq_fwd}->{name_desc} ."\t".
			$ann_enst2geneID->{$refseq_fwd}->{description} ."\t";
	}

	else {
		print NDG "NA\tNA\tNA\tNA\tNA\t";
	}


	### Print rev gene annotation
	print NDG "$refseq_rev\t$refseq_rev_dist_txStart\t$refseq_rev_dist_cdsStart\t". 
			$ann_UCSC_RefSeqGenes_rev->{$refseq_rev}->{name2} ."\t";

	if( ($refseq_rev ne "NA") && (exists $ann_NCBI_gene2ensembl->{$refseq_rev}->{GeneID}) ) {
		my $geneID = $ann_NCBI_gene2ensembl->{$refseq_rev}->{GeneID};

		print NDG $geneID ."\t".
			$ann_NCBI_gene2ensembl->{$refseq_rev}->{Ensembl_rna_id} ."\t".
			$ann_NCBI_gene2ensembl->{$refseq_rev}->{Ensembl_gene_id} ."\t".
			$ann_NCBI_gene_info->{$geneID}->{Synonyms} ."\t".
			$ann_NCBI_gene_info->{$geneID}->{description} ."\n";
	}
	elsif( ($refseq_rev ne "NA") && (exists $ann_enst2geneID->{$refseq_rev}->{GeneID}) ) {
		my $geneID = $ann_enst2geneID->{$refseq_rev}->{GeneID};

		print NDG $geneID ."\t".
			$ann_enst2geneID->{$refseq_rev}->{enst_name} ."\t".
			$ann_enst2geneID->{$refseq_rev}->{ensg_name} ."\t".
			$ann_enst2geneID->{$refseq_rev}->{name_desc} ."\t".
			$ann_enst2geneID->{$refseq_rev}->{description} ."\n";
	}
	else {
		print NDG "NA\tNA\tNA\tNA\tNA\n";
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


