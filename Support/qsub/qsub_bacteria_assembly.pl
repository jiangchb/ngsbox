#!/usr/bin/perl
use warnings;
#use strict;

##Read sequences from input folder and define options

my $usage = "\n$0 infolder outfolder qsubName num_threads\n\n";
my $infolder  = shift or die $usage;
my $outfolder = shift or die $usage;
my $qsub_name = shift or die $usage;
my $t   = shift or die $usage;
my $PWD=`pwd`;
chomp $PWD;

my @files = glob($infolder . "/*1.fastq");


foreach my $read1 (@files) {

	# read 2 filename
	my $read2 = $read1;
	$read2 =~ s/1.fastq/2.fastq/g;
	
	my @filepath = split("/", $read1);
	my @filepath2 = split("/", $read2);
	my $fileleaf1 = $filepath[$#filepath];
	my $fileleaf2 = $filepath2[$#filepath2];
	my $name = substr($fileleaf1, 0, 9);

if(! -e "$outfolder/$name") {
		mkdir "$outfolder/$name" or die "Cannot create output directory $outfolder/$name";
		mkdir "$outfolder/$name/assembly" or die "Cannot create output directory $outfolder/$name/assembly";
		mkdir "$outfolder/$name/alignment" or die "Cannot create output directory $outfolder/$name/alignment";
	}
	else {
		print STDERR "Runfolder $outfolder/$name already exists. Will only update qsub file\n";
	}

	open OUT, ">$outfolder/$name/$qsub_name" or die "Cannot create qsub file $outfolder/$name/$qsub_name";

my @qsub = ("#!/bin/bash

# \$ -e $outfolder/$name/
# \$ -o $outfolder/$name/

NAME=$name
READ1=$read1
READ2=$read2
OUTF=$PWD/$outfolder/$name
OUTA=$PWD/$outfolder/$name/assembly
OUTS=$PWD/$outfolder/$name/scaffolding
OUTD=$PWD/$outfolder/$name/alignment
REF2=/users/GD/resource/Thermotoga/NC_000853.SB/NC_000853.SB.fasta
REF1=/users/GD/resource/Pseudomonas/Pseudomonas_aeruginosa_PA7_uid58627/NC_009656.fasta
REF3=/users/GD/resource/Pseudomonas/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fasta
REF4=/users/GD/resource/Pseudomonas/Pseudomonas_aeruginosa_UCBPP_PA14_uid57977/NC_008463.fasta
PERL=/soft/bin/perl
SGA=/soft/bin/sga
SSPACE=/users/GD/projects/Plants/Ath/tools/ngopt-read-only/bin/SSPACE/SSPACE
GAPFILLER=/users/so/lzapata/tools/GapFiller_v1-10_linux-x86_64/GapFiller.pl
MUMMER=/users/GD/tools/mummer/MUMmer3.23
GAGE=/users/GD/tools/gage/getCorrectnessStats2.sh
SOAP=/users/GD/tools/soap/SOAPdenovo-V1.05/SOAPdenovo31mer
MV=/bin/mv
MKDIR=/bin/mkdir
THREAD=$t
PWD=$PWD
##
##Filter reads with illumina Bad quality Flag

\$PERL /users/so/lzapata/newscripts/Filter/filterbyIlluminaflag.pl \$READ1 > \$OUTF/\$NAME.r1
\$PERL /users/so/lzapata/newscripts/Filter/filterbyIlluminaflag.pl \$READ2 > \$OUTF/\$NAME.r2

##
##Call preprocessing pipeline using SGA, (check for adapter used to construct the library and place it under your home_directory)
\$PERL /users/so/lzapata/newscripts/Correction/correct_paired_end_bacteria.pl \$OUTF/\$NAME.r1 \$OUTF/\$NAME.r2 \$NAME \$THREAD \$OUTF/

\$MKDIR \$OUTA
\$MKDIR \$OUTD
\$MKDIR \$OUTS

\$MV \$OUTF/\$NAME.r1.pp.ec.fastq \$OUTF/\$NAME.finalr1.fastq
\$MV \$OUTF/\$NAME.r2.pp.ec.fastq \$OUTF/\$NAME.finalr2.fastq

##
##Print config File for SOAP to run
echo '#maximal read length \
max_rd_len=100 \
[LIB] \
avg_ins=340 \
reverse_seq=0 \
asm_flags=3 \
rd_len_cutoff=100 \
rank=1 \
pair_num_cutoff=3 \
map_len=25 \
q1='\$OUTF/\$NAME'.finalr1.fastq\
q2='\$OUTF/\$NAME'.finalr2.fastq' > \$OUTA/\$NAME.config
##


##Assemble sequences using SOAP
\$SOAP all -s \$OUTA/\$NAME.config -K 31 -p \$THREAD -R -o \$OUTA/\$NAME.soap '&>' \$OUTA/\$NAME.soap.log

##
###Print config file for SSPACE
echo 'Lib1 '\$OUTF/\$NAME'.finalr1.fastq '\$OUTF/\$NAME'.finalr2.fastq 340 0.3 0' > \$OUTS/\$NAME.sspace.lib 
##

##Do Rescaffolding with SSPACE
\$PERL \$SSPACE -l \$OUTS/\$NAME.sspace.lib -s \$OUTA/\$NAME.soap.scafSeq -x 0 -m 32 -o 20 -t 0 -k 5 -n 15 -v 0 -d \$OUTS -b \$NAME.sspace '&>' sspace.log 

##Print Config file for GapFiller
#
echo 'Lib1GF bowtie '\$OUTF/\$NAME'.finalr1.fastq '\$OUTF/\$NAME'.finalr2.fastq 340 0.75 FR' > \$OUTS/\$NAME.gapfill.lib
#

\$PERL \$GAPFILLER -l \$OUTS/\$NAME.gapfill.lib -s \$OUTS/\$NAME.sspace.final.scaffolds.fasta -m 72 -o 2 -r 0.7 -n 6 -t 0 -i 15 -d 50 -T 4 -b \$OUTF/\$NAME.gapfilled

\$MV \$OUTF/\$NAME.gapfilled/\$NAME.gapfilled.gapfilled.final.fa \$OUTF/\$NAME.gapfilled.gapfilled.final.fa

##
##Whole Genome Alignment to reference
\$MUMMER/nucmer --prefix=\$OUTD/\$NAME_ref1 \$REF1 \$OUTF/\$NAME.gapfilled.gapfilled.final.fa
##\$MUMMER/nucmer --prefix=\$OUTD/\$NAME_ref2 \$REF2 \$OUTA/\$NAME.soap.scaffSeq
##\$MUMMER/nucmer --prefix=\$OUTD/\$NAME_ref3 \$REF3 \$OUTA/\$NAME.soap.scaffSeq
##\$MUMMER/nucmer --prefix=\$OUTD/\$NAME_ref4 \$REF4 \$OUTA/\$NAME.soap.scaffSeq
##\$MUMMER/delta-filter -o 95 -i 95 \$OUTD/\$NAME_ref1.delta > \$OUTD/\$NAME_ref1.fdelta
##\$MUMMER/show-coords -lrcT \$OUTD/\$NAME_ref1.fdelta | sort -k13 -k1n -k2n > \$OUTD/\$NAME_ref1.coords
##\$MUMMER/show-tiling -c -l 1 -i 0 -V 0 \$OUTD/\$NAME_ref1.fdelta > \$OUTD/\$NAME_ref1.tiling

##
##Genome alignment with correction GAGE
##\$GAGE \$REF1 \$OUTA/\$NAME.soap.contig \$OUTA/\$NAME.soap.scafSeq \$OUTD '&>' \$OUTD/gage.ref1.report
##\$GAGE \$REF2 \$OUTA/\$NAME.soap.contig \$OUTA/\$NAME.soap.scaffSeq
##\$GAGE \$REF3 \$OUTA/\$NAME.soap.contig \$OUTA/\$NAME.soap.scaffSeq
##\$GAGE \$REF4 \$OUTA/\$NAME.soap.contig \$OUTA/\$NAME.soap.scaffSeq

\n");
	print OUT @qsub;

	close OUT;
}
