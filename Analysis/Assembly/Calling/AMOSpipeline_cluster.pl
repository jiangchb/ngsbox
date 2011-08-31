#! /usr/bin/perl
use strict;
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


