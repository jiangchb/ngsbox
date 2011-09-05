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
#  Module: Annotation::Variants::WGA_annotate_polymorphic_genes.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

### User params
my $usage = "\n\n$0 ecotype geneGFF genesFolder partialGenesFolder centromere\n\n";
my $ecotype            = shift or die $usage;
my $geneGFF            = shift or die $usage;
my $genesFolder        = shift or die $usage;
my $partialGenesFolder = shift or die $usage;
my $centro             = shift or die $usage;

my $dbh;
&connect_to_db();

my %AGI = ();
my %indel_length_dist = ();
my %hdr_length_dist = ();

### Read centromere regions file
my %centro_regions = ();
open CENTRO, $centro or die "Cannot open $centro\n\n";
while(<CENTRO>) {
	chomp;
	my @e = split("\t", $_);
	for(my $i = $e[1]; $i <= $e[2]; $i++) {
		$centro_regions{$i} = $e[0];
	}
}


### Parse gene annotation in GFF format
open GENE, $geneGFF or die "Cannot open $geneGFF\n";
while(<GENE>) {
	chomp;
	my @e = split("\t", $_);

	if( $e[2] eq "gene" ) {

		my @ids = split(";", $e[8]);
		my $gene_id = substr($ids[0], 3);
		my $gene_chr = substr($e[0], 3, 1);

		if( $gene_id =~ /AT\dG.*/ ) {
			if( (! exists $centro_regions{$e[3]}) || ($centro_regions{$e[3]} != $gene_chr) ) {
				my %AGI_entry = (AAchange => 0, newstop => 0, loststop => 0, indel_pos => -1, FS => 0, nonFS => 0,  
						microindel => 0, microdel => 0, microins => 0, hdr_genes => 0, hdr_cds => 0, full => 0, partial => 0);
				$AGI{$gene_id} = \%AGI_entry;
			}
		}
	}
}


### Read full gene alignments
my @geneAlis = glob("$genesFolder/AT*G*");

foreach my $geneAli ( @geneAlis ) {

	my @subpath = split("/", $geneAli);
	my $subleaf = $subpath[$#subpath];
	my $gene_id = substr($subleaf, 0, 9);
	if( exists $AGI{$gene_id}) {
		$AGI{$gene_id}{full} = 1;
	}	
}


### Read partial gene alignments
my @partialGeneAlis = glob("$partialGenesFolder/AT*G*");

foreach my $partialGeneAli ( @partialGeneAlis ) {

	my @subpath = split("/", $partialGeneAli);
	my $subleaf = $subpath[$#subpath];
	my $gene_id = substr($subleaf, 0, 9);
	if( exists $AGI{$gene_id}) {
		$AGI{$gene_id}{partial} = 1;
	}	
}


### Load SNP
#ecotype | chromosome | position | reference | variant | seq_type   | AGI  | isoform | CDSpos | codonpos | syn  | refAA | varAA
my $q ="SELECT 	AGI, refAA, varAA 
	FROM 	assembly_ann_snp
	WHERE 	ecotype = '$ecotype' &&
		seq_type = 'CDS'
	ORDER BY chromosome, position
";
my $sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {

	if( exists $AGI{$ref->{AGI}}) {
		if($ref->{varAA} eq "*") {
			$AGI{$ref->{AGI}}{newstop} = 1;
		}
		elsif($ref->{refAA} eq "*") { 
			$AGI{$ref->{AGI}}{loststop} = 1;
		}
		else {
			$AGI{$ref->{AGI}}{AAchange} = 1;
		}
	}
}


### Load micro-deletions
#ecotype | chromosome | begin | end | seq | seq_type | AGI
$q = "	SELECT 	AGI, seq
	FROM	assembly_ann_del
	WHERE 	ecotype = '$ecotype' &&
		seq_type = 'CDS'
	ORDER by chromosome, begin
";
$sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {
	if(exists $AGI{$ref->{AGI}}) {
		my $len = length($ref->{seq});
		$indel_length_dist{$len}++;
		
		$AGI{$ref->{AGI}}{microdel} = 1;
		$AGI{$ref->{AGI}}{microindel} += $len;

		if($len % 3 == 0) {
			$AGI{$ref->{AGI}}{nonFS}++;
		}
		else {
			$AGI{$ref->{AGI}}{FS}++;
		}
	}
}


### Load micro-insertions
#ecotype | chromosome | begin | end | seq | seq_type | AGI
$q = "  SELECT  *
	FROM    assembly_ann_ins
	WHERE   ecotype = '$ecotype' &&
		seq_type = 'CDS'
	ORDER by chromosome, begin
";
$sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {
	if(exists $AGI{$ref->{AGI}}) {
		my $len = length($ref->{seq});
		$indel_length_dist{$len}++;
		
		$AGI{$ref->{AGI}}{microins} = 1;
		$AGI{$ref->{AGI}}{microindel} -= $len;

		if($len % 3 == 0) {
			$AGI{$ref->{AGI}}{nonFS}++;
		}
		else {
			$AGI{$ref->{AGI}}{FS}++;
		}	
	}
}


### Load HDRs
#ecotype | chromosome | begin | end | indel_type | seq_type | gene_absence | utr_absence | cds_absence
$q = "  SELECT  chromosome, begin, end, indel_type, seq_type, gene_absence, cds_absence
	FROM    assembly_ann_hdr
	WHERE   ecotype = '$ecotype' &&
		seq_type != 'intergenic'
	ORDER by chromosome, begin
";
$sth = $dbh->prepare($q);
$sth->execute();

while(my $ref = $sth->fetchrow_hashref()) {
	my $len = $ref->{end} - $ref->{begin} + 1;
	$hdr_length_dist{$len}++;

	my @gene_absence = split(",", $ref->{gene_absence});
	foreach my $gene ( @gene_absence ) {
		if(exists $AGI{$gene}) {
			$AGI{$gene}{hdr_genes} = 1;
		}
	}
	
	my @cds_absence = split(",", $ref->{cds_absence});
	foreach my $cds ( @cds_absence ) {
		if(exists $AGI{$cds}) {
			$AGI{$cds}{hdr_cds} = 1;
		}
	}
}


### Final statistics
# AAchange => 0, newstop => 0, loststop => 0, microindel => 0, microdel => 0, microins => 0, hdr => 0, full => 0, partial => 0
my $total_genes        = 0;
my $accessible         = 0;
my $accessible_full    = 0;
my $accessible_partial = 0;
my $hdr_gene           = 0;
my $hdr_cds            = 0;
my $major              = 0;
my $minor              = 0;
my $conserved          = 0;
my $majordel           = 0;
my $majorins           = 0;
my $minordel           = 0;
my $minorins           = 0;
my $newstop            = 0;
my $loststop           = 0;
my $aachange           = 0;
my $compensated        = 0;
my $nonFS              = 0;

my %gene_length_change_dist = ();

foreach my $gene (keys %AGI) {

	# Accessibility
	$total_genes++;

	if($AGI{$gene}{full} == 1) {
		$accessible++;
		$accessible_full++;
	}
	elsif($AGI{$gene}{partial} == 1) {
		$accessible++;
		$accessible_partial++;
	}

	# Absence/Presence polymorphism or strongly differing alleles
	if($AGI{$gene}{hdr_genes} == 1) {
		$hdr_gene++;
	}
	if($AGI{$gene}{hdr_cds} == 1) {
		$hdr_cds++;
	}

	# Genes with full alignment
	if($AGI{$gene}{full} == 1) {

		# Potential major change -> lost or new stop or indel not mod(3) = 0
		if( ($AGI{$gene}{newstop} == 1) || ($AGI{$gene}{loststop} == 1) || ( ($AGI{$gene}{microindel} % 3) != 0) ) {
			
			$major++;
			
			if($AGI{$gene}{microdel} == 1) { $majordel++; }
			if($AGI{$gene}{microins} == 1) { $majorins++; }
			if($AGI{$gene}{newstop} ==  1) { $newstop++;  }
			if($AGI{$gene}{loststop} == 1) { $loststop++; }
		}

		# Minor change -> AA change or indel mod(3) = 0
		elsif( ($AGI{$gene}{AAchange} == 1) || ($AGI{$gene}{microdel} == 1) || ($AGI{$gene}{microins} == 1) ) {
			
			$minor++;

			if($AGI{$gene}{microdel} == 1) { $minordel++; }
			if($AGI{$gene}{microins} == 1) { $minorins++; }
			if($AGI{$gene}{AAchange} == 1) { $aachange++; }

			# Compensatory indels
			if( ($AGI{$gene}{microdel} == 1) || ($AGI{$gene}{microins} == 1) ) {
				if($AGI{$gene}{FS} > 1) {
					$compensated++;
				}
				elsif($AGI{$gene}{nonFS} > 0) {
					$nonFS++;
				}
				else {
					print STDERR "PROBLEMS\n";
				}
			}
		}

		# Conserved
		else {
			$conserved++;
		}

		# Gene length change distribution
		if( ($AGI{$gene}{microdel} == 1) || ($AGI{$gene}{microins} == 1) ) { 
			$gene_length_change_dist{ abs($AGI{$gene}{microindel}) }++;
		}
	}
}

print "

Total genes: $total_genes\n
Accessible: $accessible\n
Accessible (full alignment): $accessible_full\n
Accessible (partial alignment): $accessible_partial\n
Absence/Presence polymorphism in genes: $hdr_gene\n
Absence/Presence polymorphism in CDS: $hdr_cds\n\n
Major change: $major\n
\tmajor del: $majordel\n
\tmajor ins: $majorins\n
\tnew stop: $newstop\n
\tlost stop: $loststop\n\n
Minor change: $minor\n
\tAA change: $aachange\n
\tminor del: $minordel\n
\tminor ins: $minorins\n
\tcompensatory indels: $compensated\n
\tin frame indels: $nonFS\n
Conserved:$conserved\n
\n";


open INDDIST, ">indel_length_dist.txt" or die "Cannot open output file\n";
foreach my $len (sort {$a<=>$b} keys %indel_length_dist) {
	print INDDIST "$len\t" . $indel_length_dist{$len} . "\n";
}
close INDDIST;

open HDRDIST, ">hdr_length_dist.txt" or die "Cannot open output file\n";
foreach my $len (sort {$a<=>$b} keys %hdr_length_dist) {
	print HDRDIST "$len\t" . $hdr_length_dist{$len} . "\n";
}
close HDRDIST;

open GENEINDDIST, ">gene_length_change_dist.txt" or die "Cannot open output file\n";
foreach my $len (sort {$a<=>$b} keys %gene_length_change_dist) {
	print GENEINDDIST "$len\t" . $gene_length_change_dist{$len} . "\n";
}


exit(0);


#####################################################
### Connects to a database and returns databaseHandle
sub connect_to_db
{
        my $databaseName = "ath_pe";
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

