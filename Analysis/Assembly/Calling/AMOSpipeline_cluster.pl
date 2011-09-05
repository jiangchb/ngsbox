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
#  Module: Analysis::Assembly::Calling::AMOSpipeline_cluster.pl
#  Purpose:
#  In:
#  Out:
#

use Cwd;

my $usage = "$0 shorebinary amosbinary AMOS_T minContigSize ScaffoldFolder [ScaffoldFolder [ScaffoldFolder [...]]]\n";
my $shore_bin = shift or die $usage; 
my $amos_bin = shift or die $usage; 
my $amos_T = shift or die $usage; 
my $contig_size = shift or die $usage; 

#my $a_param = "velvet,superlocas,euler,abyss";
#my $a_param = "superlocas,euler,abyss";
my $a_param = "velvet";
#my $a_param = "velvet,euler,abyss";
#my $a_param = "velvet,superlocas,abyss";
#my $a_param = "velvet,superlocas,euler";
#my $a_param = "velvet,superlocas";

my @folders = @ARGV;

my $basedir = getcwd;

foreach my $folder (@folders) {
	chdir($basedir);
	chdir($folder."/AMOScmp_batches");
	my @batches = glob("AMOS_batch_*");

	foreach my $batch (@batches) {
		chdir($batch);
		if (not -e "contigs.fasta" and $batch !~ m/test/) {
			print $basedir."/".$folder."/AMOScmp_batches/".$batch."\n";
			system("rm -rf contigs.afg contigs.bnk contigs.cluster contigs.conflict contigs.delta contigs.layout contigs.runAmos.log contigs.seq left_over_contigs contigs.ntref contigs.mgaps nucmer.error run_amos.sh*");
			system("$shore_bin convert -a $a_param -c Contig2AFG -s 5000 -m $contig_size -i contigs.fa -o contigs.afg");
		
			# write simple cluster script
			my $jobdir = getcwd;
			open CS, ">run_amos.".$amos_T.".sh";
			print CS "#!/bin/sh\n";
			print CS "JOBDIR=".$jobdir."\n";
			print CS "JOBOUT=\"\$JOBDIR\"\n";
			print CS "if [ ! -e contigs.fasta ]; then $amos_bin contigs -D REF=ref.fasta > \$JOBOUT/amos.".$amos_T.".out; fi\n";
			close CS;

			system("qsub -R y -p 1000 -cwd -l h_vmem=90G run_amos.".$amos_T.".sh")
		}
		chdir("..");
	}
}


