#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

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
#  Module: Support::qsub::qsub_batch_exome_fastq_SE.pl
#  Purpose:
#  In:
#  Out:
#
sub usage { print "\n$0 \n usage:\n",
	   "--infolder \t folder containing the raw sequence data in fastq format \n",
	   "--outfolder \t folder to write the analysis to \n",
	   "--qsubname \t name of the script you want to start lateron \n",
	   "--max_coverage \t used in SNP filtering with samtools (400 works well) \n",
	   "--namestart \t  start of a substring in the read file name (first letter is numbered 1) \n",
	   "--namelength \t length of the part of the filenames which should be taken as sample (and folder) name \n",
	   "--readextension \t describes the name of the read files (e.g. fastq.gz)\n",
           "--enrichment \t which enrichment kit was used? (Agilent: a35s - 35MB standard, a35e - 35MB extended, a50 - 50MB standard, Nimblegene: n50 - 50MB version 3)\n",
	   "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 8]\n",
	   "--help \t\t show help \n";
	   
	   }


my $infolder;
my $outfolder;
my $qsub_name;
my $max_cov; 
my $nameStart = 'NA';
my $nameLength;
my $readextension;
my $enrichment;
my $cpu = 8;
my $help = 0;


my %enrichmentkits = ('a35s' => '/users/GD/resource/human/probesets/agilent/35MB_standard/shore_format',
		      'a35e' => '/users/GD/resource/human/probesets/agilent/35MB_extended/shore_format',
		      'a50' => '/users/GD/resource/human/probesets/agilent/50MB/shore_format',
		      'n50' => '/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/shore_format/');


GetOptions("infolder=s" => \$infolder, "outfolder=s" => \$outfolder, "qsubname=s" => \$qsub_name, "max_coverage=i" => \$max_cov, "nameStart=s" => \$nameStart, "nameLength=s" => \$nameLength, "readextension=s" => \$readextension, "cpu=i" => \$cpu, "enrichment=s" => \$enrichment, "help=s" => \$help);


unless($infolder && $outfolder && $qsub_name && $max_cov && $nameStart ne 'NA' && $nameLength && $readextension && $enrichment && $help == 0) {
	usage;
	exit;
}

unless(defined($enrichmentkits{$enrichment})) {
	print "\n###########################\n enrichment method unknown \n###########################\n";
	usage;
	exit;
}

my @files = glob($infolder . "/*$readextension");

foreach my $file (@files) {

	my @filepath = split("/", $file);
	my $fileleaf = $filepath[$#filepath];


	my $name = substr($fileleaf, $nameStart - 1, $nameLength);

	if(! -e "$outfolder/$name") {
		mkdir "$outfolder/$name" or die "Cannot create output directory $outfolder/$name";
	}
	else {
		print STDERR "Runfolder $outfolder/$name already exists. Will only update qsub file\n"
	}

	open OUT, ">$outfolder/$name/$qsub_name" or die "Cannot create qsub file $outfolder/$name/$qsub_name";

my $qsubNname; 
if ($name =~ /^\d/) {
	$qsubNname = 's'.$name;
}
else {
	$qsubNname = $name;
}# qsub doesn't like jobs starting with a digit

my @qsub = ("#!/bin/bash

#\$ -N $qsubNname
#\$ -e $outfolder/$name/
#\$ -o $outfolder/$name/


source /users/so/sossowski/.bashrc
export TMPDIR=/users/GD/projects/HumanDisease/tmp
export _JAVA_OPTIONS=-Djava.io.tmpdir=/users/GD/projects/HumanDisease/tmp
export PATH=/users/GD/tools/annovar/annovar_2011Nov20/:\$PATH


NAME=$name
FASTQ=$file
OUTF=$outfolder/$name
REF=/users/GD/resource/human/hg19/bwa/hg19.fasta
SHOREREF=/users/GD/resource/human/hg19/shore/hg19.fasta.shore
EXOME=$enrichmentkits{$enrichment}
KNOWNINDEL=/users/GD/resource/human/hg19/databases/dbSNP/dbIndel132_20101103.vcf
KNOWN=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp135_20111104.vcf
BWA=/users/GD/tools/bwa/bwa-0.5.10/bwa
GATK=/users/GD/tools/GATK/GATK_src_1.4-15-gcd43f01/dist/GenomeAnalysisTK.jar
SAMTOOLS=/soft/bin
ANNOVAR=/users/GD/tools/annovar/annovar_2011Nov20
SHORE=/users/GD/tools/shore/shore
NGSBOX=/users/GD/tools/ngsbox
RSCRIPT=/soft/bin/Rscript
 

### Align reads with bwa
 \$BWA aln -k 2 -i 5 -q -1 -t $cpu -R 0 -n 6 -o 1 -e 20 -l 28 -f \$OUTF/\$NAME.sai \$REF \$FASTQ



### Create SAM files
if [ -s \$OUTF/\$NAME.sai ];
then
   echo BWA samse
   \$BWA samse -r \"\@RG\\tID:\$NAME\\tSM:\$NAME\" -n 10 -f \$OUTF/\$NAME.sam \$REF \$OUTF/\$NAME.sai \$FASTQ
else
   echo \$NAME.sai
   exit
fi



### Convert SAM to BAM
if [ -s \$OUTF/\$NAME.sam ];
then
   echo Convert SAM to BAM
   \$SAMTOOLS/samtools view -b -S -o \$OUTF/\$NAME.bam \$OUTF/\$NAME.sam
else
   echo \$NAME.sam not found
   exit
fi



### Sort BAM file
if [ -s \$OUTF/\$NAME.bam ];
then
   echo Sort BAM
   \$SAMTOOLS/samtools sort \$OUTF/\$NAME.bam \$OUTF/\$NAME.sort
else
   echo \$NAME.bam not found
   exit
fi



### Index sorted BAM file
if [ -s \$OUTF/\$NAME.sort.bam ];
then
   echo Index sorted BAM file
   \$SAMTOOLS/samtools index \$OUTF/\$NAME.sort.bam
else
   echo \$NAME.sort.bam not found
   exit
fi



### Local Re-alignment and 'fake dm'
if [ -s \$OUTF/\$NAME.sort.bam.bai ];
then
   echo Local Re-alignment
   java -Xmx4g -jar \$GATK -T RealignerTargetCreator -R \$REF -I \$OUTF/\$NAME.sort.bam -o \$OUTF/\$NAME.intervals -known \$KNOWNINDEL --minReadsAtLocus 6 --maxIntervalSize 200
   java -Xmx4g -jar \$GATK -T IndelRealigner -R \$REF -I \$OUTF/\$NAME.sort.bam -targetIntervals \$OUTF/\$NAME.intervals -o \$OUTF/\$NAME.realigned.dm.bam -known \$KNOWNINDEL --maxReadsForRealignment 10000 --consensusDeterminationModel USE_SW
else
   echo \$NAME.sort.bam.bai not found
   exit
fi


# ---->>>> end of difference


### Base quality recalibration
if [ -s \$OUTF/\$NAME.realigned.dm.bam ];
then
   echo Base quality recalibration
   java -Xmx4g -jar \$GATK -T CountCovariates -nt $cpu --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile \$OUTF/recal_data.csv -R \$REF -I \$OUTF/\$NAME.realigned.dm.bam -knownSites \$KNOWN
   java -Xmx4g -jar \$GATK -T TableRecalibration --default_platform illumina -R \$REF -I \$OUTF/\$NAME.realigned.dm.bam -recalFile \$OUTF/recal_data.csv --out \$OUTF/\$NAME.realigned.dm.recalibrated.bam
else
   echo \$NAME.realigned.dm.bam not found
   exit
fi



### Cleanup
rm \$OUTF/\$NAME.sam
rm \$OUTF/\$NAME.bam
rm \$OUTF/\$NAME.realigned.bam
rm \$OUTF/\$NAME.realigned.bai
rm \$OUTF/\$NAME.realigned.dm.bam
rm \$OUTF/\$NAME.realigned.dm.bam.bai



### GATK: Call SNPs and Indels with the GATK Unified Genotyper
if [ -s \$OUTF/\$NAME.realigned.dm.recalibrated.bam ];
then
   echo GATK: Call SNPs and Indels with the GATK Unified Genotyper
   java -Xmx4g -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.snps.raw.vcf -glm SNP
   java -Xmx4g -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.indel.raw.vcf -glm INDEL
else
   echo \$NAME.realigned.dm.recalibrated.bam not found
   exit
fi

if [ ! -s \$OUTF/GATK.snps.raw.vcf ];
then
   echo GATK.snps.raw.vcf not found
   exit
fi



### MPILEUP: Call SNPs and Indels
 \$SAMTOOLS/samtools mpileup -uf \$REF \$OUTF/\$NAME.realigned.dm.recalibrated.bam | \$SAMTOOLS/bcftools view -bcg - > \$OUTF/MPILEUP.variant.raw.bcf
 \$SAMTOOLS/bcftools view \$OUTF/MPILEUP.variant.raw.bcf | \$SAMTOOLS/vcfutils.pl varFilter -d5 -D$max_cov -W 20 > \$OUTF/MPILEUP.variant.raw.vcf
egrep \"INDEL|#\" \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.indel.raw.vcf
grep -v INDEL \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.snps.raw.vcf

if [ ! -s \$OUTF/MPILEUP.snps.raw.vcf ];
then
   echo MPILEUP.snps.raw.vcf not found
   exit
fi



### SHORE: Prepare format map.list
mkdir \$OUTF/shore
 \$SHORE convert --sort -r \$REF -n 6 -g 1 -e 20 -s Alignment2Maplist \$OUTF/\$NAME.realigned.dm.recalibrated.bam \$OUTF/shore/map.list.gz

if [ ! -s \$OUTF/shore/map.list.gz ];
then
   echo shore/map.list.gz not found
   exit
fi



### SHORE: compute coverage plot in GFF format for browsers
 \$SHORE coverage -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/CoverageAnalysis



### SHORE: Compute enrichment
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus200 -f \$EXOME/Exome_Array_plus200.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus150 -f \$EXOME/Exome_Array_plus150.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus100 -f \$EXOME/Exome_Array_plus100.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus50 -f \$EXOME/Exome_Array_plus50.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus0 -f \$EXOME/Exome_Array_plus0.bed -H 1,1 -k



### SHORE: Enrichment plots
grep enriched \$OUTF/shore/Count_SureSelect_plus150/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus150/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_depleted.txt
grep enriched \$OUTF/shore/Count_SureSelect_plus150/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_count_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus150/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_count_depleted.txt
### plot data
\$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus150/

grep enriched \$OUTF/shore/Count_SureSelect_plus0/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus0/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_depleted.txt
grep enriched \$OUTF/shore/Count_SureSelect_plus0/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_count_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus0/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_count_depleted.txt
### plot data
\$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus0/

grep enriched \$OUTF/shore/Count_SureSelect_plus50/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus50/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_depleted.txt
grep enriched \$OUTF/shore/Count_SureSelect_plus50/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_count_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus50/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_count_depleted.txt
### plot data
\$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus50/

grep enriched \$OUTF/shore/Count_SureSelect_plus100/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus100/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_depleted.txt
grep enriched \$OUTF/shore/Count_SureSelect_plus100/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_count_enriched.txt
grep depleted \$OUTF/shore/Count_SureSelect_plus100/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_count_depleted.txt
### plot data
\$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus100/



### SHORE: Call SNPs and Indels
 \$SHORE qVar -n \$NAME -f \$SHOREREF -o \$OUTF/shore/Variants -i \$OUTF/shore/map.list.gz -s /users/GD/tools/shore/scoring_matrices/scoring_matrix_het.txt -E \$OUTF/shore/Count_SureSelect_plus150/meancov.txt -e -c 4 -d 4 -C $max_cov -r 3 -q 10 -Q 15 -a 0.25 -b 6 -y -v

if [ ! -s \$OUTF/shore/Variants/ConsensusAnalysis ];
then
   echo shore/Variants/ConsensusAnalysis is empty
   exit
fi



### Clean up
rm -r \$OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
gzip -9 \$OUTF/shore/Variants/ConsensusAnalysis/reference.shore




### Filter and compare SNP calls from 3 different pipelines
# Filtering
mkdir \$OUTF/SNP_Intersection

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf --variant \$OUTF/GATK.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 25.0 || QD < 4.0 || HRun > 9\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.vcf --variant \$OUTF/MPILEUP.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 15.0 || DP < 5 || DP > $max_cov\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter.pl \$OUTF/shore/Variants/ConsensusAnalysis/snp.vcf $max_cov > \$OUTF/SHORE.snps.raw.vcf

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/SHORE.snps.filtered.vcf --variant \$OUTF/SHORE.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"QUAL < 20.0 || DP < 5 || DP > $max_cov\" --filterName CRG

# greping
grep -v \"CRG\" \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
grep -v \"CRG\" \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf
grep -v \"CRG\" \$OUTF/SNP_Intersection/SHORE.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf

# Correct sample names in VFC files
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-MPILEUP/\" \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-SHORE/\" \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf

if [[ ! ( -s \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf && -s \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf && -s \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf ) ]];
then
   echo GATK.snps.filtered.cleaned.vcf or MPILEUP.snps.filtered.cleaned.vcf or SHORE.snps.filtered.cleaned.vcf not found
   exit
fi

# Intersecting
java -jar -Xmx4g \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:SHORE \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf -V:GATK \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -V:MPILEUP \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf -priority GATK,MPILEUP,SHORE -o \$OUTF/SNP_Intersection/merged.vcf

# Evaluation
java -jar -Xmx4g \$GATK -T VariantEval -R \$REF --dbsnp \$KNOWN -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/SNP_Intersection/report.all.txt --eval \$OUTF/SNP_Intersection/merged.vcf -l INFO

# Annotate Enrichment
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$EXOME/Exome_Array_plus150.bed \$OUTF/SNP_Intersection/merged.vcf > \$OUTF/SNP_Intersection/merged.all.vcf

# Evaluate calls on enriched regions
grep -v \"NOTENRICHED\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.enriched.vcf

java -jar -Xmx4g \$GATK -T VariantEval -R \$REF --dbsnp \$KNOWN -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/SNP_Intersection/report.enriched.txt --eval \$OUTF/SNP_Intersection/merged.enriched.vcf -l INFO




### Filter and compare indel calls from 3 different pipelines
# Filtering
mkdir \$OUTF/Indel_Intersection

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant \$OUTF/GATK.indel.raw.vcf --filterExpression \"MQ < 30.0 || QUAL < 20.0 || MQ0 > 5 || QD < 4.0 || HRun > 9\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.vcf --variant \$OUTF/MPILEUP.indel.raw.vcf --filterExpression \"MQ < 30.0 || QUAL < 10.0 || DP < 5 || DP > $max_cov\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/SHORE.indel.filtered.vcf --variant \$OUTF/shore/Variants/ConsensusAnalysis/indels.vcf --filterExpression \"QUAL < 2.0 || DP < 4 || DP > $max_cov || RE > 1.3\" --filterName CRG

# greping
grep -v \"CRG\" \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf > \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
grep -v \"CRG\" \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.vcf > \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf
grep  -v \"CRG\" \$OUTF/Indel_Intersection/SHORE.indel.filtered.vcf | grep -v \"SHOREFILTER\" > \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf

# Correct sample names in VFC files
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-MPILEUP/\" \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-SHORE/\" \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf

if [[ ! ( -s \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf && -s \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf && -s \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf ) ]];
then
   echo GATK.indel.filtered.cleaned.vcf or MPILEUP.indel.filtered.cleaned.vcf or SHORE.indel.filtered.cleaned.vcf not found
   exit
fi
 
# Intersecting
java -jar -Xmx4g \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:SHORE \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf -V:GATK \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -V:MPILEUP \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf -priority GATK,MPILEUP,SHORE -o \$OUTF/Indel_Intersection/merged.vcf

# Evaluation
java -jar -Xmx4g \$GATK -T VariantEval -R \$REF --dbsnp \$KNOWN -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/Indel_Intersection/report.all.txt --eval \$OUTF/Indel_Intersection/merged.vcf -l INFO

# Annotate Enrichment
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$EXOME/Exome_Array_plus150.bed \$OUTF/Indel_Intersection/merged.vcf > \$OUTF/Indel_Intersection/merged.all.vcf

# Evaluate calls on enriched regions 
grep -v \"NOTENRICHED\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.enriched.vcf

java -jar -Xmx4g \$GATK -T VariantEval -R \$REF --dbsnp \$KNOWN -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/Indel_Intersection/report.enriched.txt --eval \$OUTF/Indel_Intersection/merged.enriched.vcf -l INFO




### Annotate SNPs with ANNOVAR: Intersection, all three tools predict SNP
mkdir \$OUTF/SNP_Intersection/AnnovarIntersection
egrep \"Intersection|#\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.intersection.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.intersection.vcf > \$OUTF/SNP_Intersection/AnnovarIntersection/snps.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/SNP_Intersection/AnnovarIntersection/sum \$OUTF/SNP_Intersection/AnnovarIntersection/snps.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/


### Annotate SNPs with ANNOVAR: Partial union, at least 2 tools predict SNP
mkdir \$OUTF/SNP_Intersection/AnnovarPartialUnion
egrep \"Intersection|GATK-SHORE|MPILEUP-SHORE|GATK-MPILEUP|#\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.union.vcf > \$OUTF/SNP_Intersection/AnnovarPartialUnion/snps.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/SNP_Intersection/AnnovarPartialUnion/sum \$OUTF/SNP_Intersection/AnnovarPartialUnion/snps.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/


### Annotate SNPs with ANNOVAR: Union, any tool predicts SNP
mkdir \$OUTF/SNP_Intersection/AnnovarUnion
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/AnnovarUnion/snps.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/SNP_Intersection/AnnovarUnion/sum \$OUTF/SNP_Intersection/AnnovarUnion/snps.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/



### Annotate Indels with ANNOVAR: Intersection
mkdir \$OUTF/Indel_Intersection/AnnovarIntersection
egrep \"Intersection|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.union.vcf > \$OUTF/Indel_Intersection/AnnovarIntersection/indels.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/Indel_Intersection/AnnovarIntersection/sum \$OUTF/Indel_Intersection/AnnovarIntersection/indels.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/


### Annotate Indels with ANNOVAR: Partial Union, must include SHORE or at least two tools
mkdir \$OUTF/Indel_Intersection/AnnovarUnionShore
# egrep \"Intersection|GATK-SHORE|MPILEUP-SHORE|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
egrep \"Intersection|GATK-MPILEUP|SHORE|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.union.vcf > \$OUTF/Indel_Intersection/AnnovarUnionShore/indels.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/Indel_Intersection/AnnovarUnionShore/sum \$OUTF/Indel_Intersection/AnnovarUnionShore/indels.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/


### Annotate Indels with ANNOVAR: Union, indels predicted by any tool
mkdir \$OUTF/Indel_Intersection/AnnovarUnion
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/AnnovarUnion/indels.ann
 \$ANNOVAR/custom_summarize_annovar.pl --buildver hg19 --outfile \$OUTF/Indel_Intersection/AnnovarUnion/sum \$OUTF/Indel_Intersection/AnnovarUnion/indels.ann --ver1000g 1000g2011may \$ANNOVAR/hg19/



### Clean up
rm \$OUTF/MPILEUP.variant.raw.bcf

\n");



	print OUT @qsub;

	close OUT;
}
