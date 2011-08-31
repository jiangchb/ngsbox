#! /usr/bin/perl
use strict;
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


