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
#  Module: Analysis::Assembly::Calling::AMOSpipeline_sto.pl
#  Purpose:
#  In:
#  Out:
#

use Cwd;

my $usage = "$0 shorebinary ScaffoldFolder [ScaffoldFolder [ScaffoldFolder [...]]]\n";
my $shore_bin = shift or die $usage; 

my $a_param = "velvet,abyss,euler,superlocas";

my @folders = @ARGV;

my $basedir = getcwd;

foreach my $folder (@folders) {
	chdir($basedir);
	chdir($folder."/AMOScmp_batches_test");
	my @batches = glob("AMOS_batch_*");

	foreach my $batch (@batches) {
		chdir($batch);
		if (not -e "contigs.fasta") {
			system("rm -rf amos.out contigs.afg contigs.bnk contigs.cluster contigs.conflict contigs.delta contigs.layout contigs.runAmos.log contigs.seq left_over_contigs contigs.ntref contigs.mgaps nucmer.error run_amos.sh*");

			system("$shore_bin convert -a ".$a_param." -s 5000 -m 100 Contig2AFG contigs.fa contigs.afg");
			
			system("AMOScmp_sto contigs -D REF=ref.fasta");
		}
		chdir("..");
	}
}


