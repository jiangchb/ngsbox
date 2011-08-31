#!/usr/bin/perl

use strict;
use warnings;


### User params
my $usage = "$0 ChrName MinContigLength AMOScmp_batchesFolder\n";

my $chr         = shift or die $usage;
my $min_length  = shift or die $usage;
my $folder      = shift or die $usage;


### Get all subfolders (AMOS batches)
my @batches = glob("$folder/AMOS_batch_*");


### Debug output
print "\nJoining AMOS contigs from batches:\n";
foreach my $batch (@batches) {
	if ($batch !~ m/test/) {
		print "$batch/contigs.fasta\n";
	}
}
print "\nInto merged contig file:\n$folder/contigs_$chr.fa\n\n";


### Open new contig layout file
open LAYOUTOUT, ">$folder/contigs_$chr.layout" or die "Cannot open $folder/contigs_$chr.layout\n";


### Parse contigs of each batch and store in hash
my $id = 0;
my %ctg_seq = ();
my %ctg_layout = ();

foreach my $batch (@batches) {
	if ($batch !~ m/test/) {
	my @batch_id = split("_", $batch);

	open CONTIGIN, "$batch/contigs.fasta" or die "Cannot open $batch/contigs.fasta\n";
	open LAYOUTIN, "$batch/contigs.layout" or die "Cannot open $batch/contigs.layout\n";

	my $current_layout_contig = <LAYOUTIN>;
	my @clc = split(" ", $current_layout_contig);

	while(<CONTIGIN>) {
		chomp($_);

		if (substr($_, 0, 1) eq ">") {
			$id++;

			# Print new layout file
			print LAYOUTOUT "C $id\t$clc[2]\t$clc[3]\t$clc[4]\n";
			while(<LAYOUTIN>) {
				if( substr($_,0,1) ne "C" ) {
					print LAYOUTOUT $_;
				}
				elsif( substr($_,0,1) eq "C" ) {
					@clc = split(" ", $_);
					last;
				}
			}
		}
		else {
			$ctg_seq{$id} .= $_;
		}
	}
	close CONTIGIN;
	}
}
close LAYOUTOUT;


### Print contigs
open CONTIGOUT, ">$folder/contigs_$chr.fa" or die "Cannot open $folder/contigs_$chr.fa\n";
foreach $id (sort {$a<=>$b} keys %ctg_seq) {
	if(length($ctg_seq{$id}) >= $min_length) {
		print CONTIGOUT ">$chr" . "_$id\n" . $ctg_seq{$id} . "\n";
	}
}
close CONTIGOUT;

system("perl ~/pgsp/Support/Fasta/calc_assem_stats.pl $folder/contigs_$chr.fa");

exit(0);
