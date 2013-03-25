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
#  Module: Support::qsub::qsub_batch_exome_fastq_PE.pl
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
	   "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n",
	   "--secondreadextension \t describes the name of the first read file (e.g. 2.fastq.gz or 2_sequence.txt.gz)\n",
           "--oligos \t specify MIPs enrichment oligos file\n",
	   "--regions \t specify MIPs enrichment regions file\n",
	   "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 8]\n",
	   "--help \t\t show help \n";
}


my $infolder;
my $outfolder;
my $qsub_name;
my $max_cov; 
my $nameStart = 'NA';
my $nameLength;
my $firstreadextension;
my $secondreadextension;
my $enrichment_oligo;
my $enrichment_region;
my $cpu = 8;
my $help = 0;


GetOptions("infolder=s" => \$infolder, "outfolder=s" => \$outfolder, "qsubname=s" => \$qsub_name, "max_coverage=i" => \$max_cov, "nameStart=s" => \$nameStart, "nameLength=s" => \$nameLength, "firstreadextension=s" => \$firstreadextension, "secondreadextension=s" => \$secondreadextension, "cpu=i" => \$cpu, "oligos=s" => \$enrichment_oligo, "regions=s" => \$enrichment_region, "help=s" => \$help);

unless($infolder && $outfolder && $qsub_name && $max_cov && $nameStart ne 'NA' && $nameLength && $firstreadextension && $secondreadextension && $enrichment_oligo && $enrichment_region && $help == 0) {
	usage;
	exit;
}


my @files = glob($infolder . "/*$firstreadextension");


foreach my $read1 (@files) {

	# read 2 filename
	my $read2 = $read1;
	$read2 =~ s/$firstreadextension/$secondreadextension/ge;
	
	my @filepath = split("/", $read1);
	my $fileleaf = $filepath[$#filepath];


	my $name = substr($fileleaf, $nameStart - 1 , $nameLength);

	unless (-e "$outfolder") {
		mkdir "$outfolder"  or die "Cannot create output directory $outfolder";
	}
	
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


export TMPDIR=/users/GD/projects/HumanDisease/tmp
export _JAVA_OPTIONS=-Djava.io.tmpdir=/users/GD/projects/HumanDisease/tmp
export PATH=/users/GD/tools/annovar/annovar_2011Nov20/:\$PATH


NAME=$name
READ1=$read1
READ2=$read2
OUTF=$outfolder/$name
REF=/users/GD/resource/human/hg19/bwa/hg19.fasta
SHOREREF=/users/GD/resource/human/hg19/shore/hg19.fasta.shore
OLIGO=$enrichment_oligo
REGION=$enrichment_region
KNOWNINDEL=/users/GD/resource/human/hg19/databases/dbSNP/dbIndel132_20101103.vcf
KNOWN=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp135_20111104.vcf
BWA=/users/GD/tools/bwa/bwa-0.5.10/bwa
GATK=/users/GD/tools/GATK/GATK_src_1.4-15-gcd43f01/dist/GenomeAnalysisTK.jar
SAMTOOLS=/soft/bin
ANNOVAR=/users/GD/tools/annovar/annovar_2011Nov20
SHORE=/users/GD/tools/shore/shore
NGSBOX=/users/GD/tools/ngsbox
RSCRIPT=/soft/bin/Rscript
FLASH=/users/GD/tools/Flash/FLASH/flash
FASTQC=/users/GD/tools/FastQC/FastQC-0.10.1/fastqc



### Preprocessing MIPs reads
perl \$NGSBOX/Parser/FASTQ/trim_length.pl 1 150 100 \$READ1 > \$OUTF/\$NAME.r1.fastq
perl \$NGSBOX/Parser/FASTQ/trim_length.pl 1 150 100 \$READ2 > \$OUTF/\$NAME.r2.fastq

\$FLASH \$OUTF/\$NAME.r1.fastq \$OUTF/\$NAME.r2.fastq -o \$OUTF/out -m 148 -M 148 -x 0.1

perl \$NGSBOX/Parser/FASTQ/trim_length.pl  1 152 152 \$OUTF/out.extendedFrags.fastq > \$OUTF/\$NAME.fastq

gzip -9 \$OUTF/\$NAME.fastq

\$FASTQC -o \$OUTF \$OUTF/\$NAME.fastq.gz

rm \$OUTF/\$NAME.r1.fastq
rm \$OUTF/\$NAME.r2.fastq
rm \$OUTF/out.extendedFrags.fastq
rm \$OUTF/out.notCombined_1.fastq
rm \$OUTF/out.notCombined_2.fastq



### Align reads with bwa
 \$BWA aln -k 2 -i 5 -q -1 -t $cpu -R 0 -n 6 -o 1 -e 20 -f \$OUTF/\$NAME.sai \$REF \$OUTF/\$NAME.fastq.gz



### Create SAM files
if [ -s \$OUTF/\$NAME.sai ];
then
   echo BWA samse
   \$BWA samse -r \"\@RG\\tID:\$NAME\\tSM:\$NAME\" -n 3 -f \$OUTF/\$NAME.sam \$REF \$OUTF/\$NAME.sai \$OUTF/\$NAME.fastq.gz
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



### Local Re-alignment
# if [ -s \$OUTF/\$NAME.sort.bam.bai ];
# then
#    echo Local Re-alignment
#    java -Xmx4g -jar \$GATK -T RealignerTargetCreator -R \$REF -I \$OUTF/\$NAME.sort.bam -o \$OUTF/\$NAME.intervals -known \$KNOWNINDEL --minReadsAtLocus 6 --maxIntervalSize 200
#    java -Xmx4g -jar \$GATK -T IndelRealigner -R \$REF -I \$OUTF/\$NAME.sort.bam -targetIntervals \$OUTF/\$NAME.intervals -o \$OUTF/\$NAME.realigned.bam -known \$KNOWNINDEL --maxReadsForRealignment 10000 --consensusDeterminationModel USE_SW
# else
#    echo \$NAME.sort.bam.bai not found
#    exit
# fi



### Cleanup
rm \$OUTF/\$NAME.sam
rm \$OUTF/\$NAME.bam



### GATK: Call SNPs and Indels with the GATK Unified Genotyper
if [ -s \$OUTF/\$NAME.sort.bam ];
then
   echo GATK: Call SNPs and Indels with the GATK Unified Genotyper
   java -Xmx4g -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.sort.bam -o \$OUTF/GATK.snps.raw.vcf -glm SNP
   java -Xmx4g -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.sort.bam -o \$OUTF/GATK.indel.raw.vcf -glm INDEL
else
   echo \$NAME.sort.bam not found
   exit
fi

if [ ! -s \$OUTF/GATK.snps.raw.vcf ];
then
   echo GATK.snps.raw.vcf not found
   exit
fi



### MPILEUP: Call SNPs and Indels
 \$SAMTOOLS/samtools mpileup -d 999999 -f \$REF \$OUTF/\$NAME.sort.bam > \$OUTF/\$NAME.pileup
 \$SAMTOOLS/samtools mpileup -uf \$REF \$OUTF/\$NAME.sort.bam | \$SAMTOOLS/bcftools view -bcg - > \$OUTF/MPILEUP.variant.raw.bcf
 \$SAMTOOLS/bcftools view \$OUTF/MPILEUP.variant.raw.bcf | \$SAMTOOLS/vcfutils.pl varFilter -d5 -D$max_cov -W 10 -1 0.0 -2 0.0 -3 0.0 -4 0.0 > \$OUTF/MPILEUP.variant.raw.vcf
egrep \"INDEL|#\" \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.indel.raw.vcf
grep -v INDEL \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.snps.raw.vcf

if [ ! -s \$OUTF/MPILEUP.snps.raw.vcf ];
then
   echo MPILEUP.snps.raw.vcf not found
   exit
fi



### SHORE: Prepare format map.list
mkdir \$OUTF/shore
 \$SHORE convert --sort -r \$REF -n 6 -g 1 -e 20 -s Alignment2Maplist \$OUTF/\$NAME.sort.bam \$OUTF/shore/map.list.gz

if [ ! -s \$OUTF/shore/map.list.gz ];
then
   echo shore/map.list.gz not found
   exit
fi



### SHORE: compute coverage plot in GFF format for browsers
 \$SHORE coverage -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/CoverageAnalysis



### SHORE: Compute enrichment
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_enrichment_oligo -f \$OLIGO -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_enrichment_region -f \$REGION -H 1,1 -k



### SHORE: Call SNPs and Indels
 \$SHORE qVar -n \$NAME -f \$SHOREREF -o \$OUTF/shore/Variants -i \$OUTF/shore/map.list.gz -s /users/GD/tools/shore/scoring_matrices/scoring_matrix_het.txt -t -c 4 -d 4 -C $max_cov -r 3 -q 20 -Q 15 -a 0.25 -b 6 -y -v

if [ ! -s \$OUTF/shore/Variants/ConsensusAnalysis ];
then
   echo shore/Variants/ConsensusAnalysis is empty
   exit
fi



### Clean up
# rm -r \$OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
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
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$REGION \$OUTF/SNP_Intersection/merged.vcf > \$OUTF/SNP_Intersection/merged.all.vcf

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
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$REGION \$OUTF/Indel_Intersection/merged.vcf > \$OUTF/Indel_Intersection/merged.all.vcf

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
