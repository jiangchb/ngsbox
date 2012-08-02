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
#  Module: Analysis::Assembly::WGA::Mummer_Parser.pl
#  Purpose:
#  In:
#  Out:
#

use List::Util qw(max min);


### User params
my $usage = "\n$0 contig_file scaffold_file maxDivBases outfolder\n\n";

my $ctg_file  = shift or die $usage;
my $scf_file  = shift or die $usage;
my $maxDivBp  = shift or die $usage;
my $outfolder = shift or die $usage;



### Create folders
mkdir "$outfolder" or die;
mkdir "$outfolder/Alignments" or die;


### Run NUCMER
#system("nucmer --mumreference -b $maxDivBp -g 90 -l 25 -c 65 -f --prefix=$outfolder/wgaMum $scf_file $ctg_file");
system("nucmer --mumreference -b $maxDivBp -g 50 -l 10 -c 65 --prefix=$outfolder/wgaMum $scf_file $ctg_file");
#system("nucmer --mumreference -b $maxDivBp -g 20 -l 40 -c 100 --prefix=$outfolder/wgaMum $scf_file $ctg_file");


### Run delta filter
#system("delta-filter -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");
system("delta-filter -i 80 -l 70 -u 80 -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");
#system("delta-filter -i 97 -l 100 -u 80 -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");


### Parse scaffolds file
open SCAFFOLDS, $scf_file or die "Cannot open $scf_file\n";

my $id = "";
my %scaffold_seqs = ();
while(<SCAFFOLDS>) {
        chomp;
	
	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
		$scaffold_seqs{$id} = "";
	}
	else {
		$scaffold_seqs{$id} .= $_;
	}
}


### Parse contigs file
open CONTIGS, $ctg_file or die "Cannot open $ctg_file\n";

$id = "";
my %CTG = ();

while(<CONTIGS>) {
	chomp;

	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		$CTG{$id} .= $_;
	}
}
close CONTIGS;


### Open result files and hashes
open HDRFILE, ">$outfolder/hdr.txt" or die  "Cannot open $outfolder/hdr.txt\n";


foreach my $scaffold (sort keys %scaffold_seqs) {

	my $SEQ = $scaffold_seqs{$scaffold};

	### Create alignments for each contig ID and parse alignment results (e.g. SNPs, Indels, HDRs, gene annotation)
	foreach $id (keys %CTG) {

		my $exit_code = system("show-aligns -r $outfolder/wgaMum.filter $scaffold $id > $outfolder/Alignments/alignments_$scaffold-$id.txt");
		if($exit_code != 0) {
			unlink("$outfolder/Alignments/alignments_$scaffold-$id.txt");
			next;
		}

		open ALIGNMENT, "$outfolder/Alignments/alignments_$scaffold-$id.txt" or die "Cannot open $outfolder/Alignments/alignments_$scaffold-$id.txt\n";


		### Init variables
		my $HSPcounter = 0;
		my %HSP = ();
		my $total_aligned = 0;
		my $strand = "NA";


		### Parse alignment file of one scaffold
		while(<ALIGNMENT>) {
			chomp;

			### New HSP found ----------------------------------------------------------------------
			if( $_ =~ /BEGIN/ ) {
				$HSPcounter++;
				my @e = split(" ", $_);
				
				# check strand consistency
				if($strand eq "NA") { 
					$strand = $e[9]; 
				}
				elsif( $strand ne $e[9] ) {
					$strand = "inconsistent";
				}

				my %current_HSP = ( 'scfbeg' => $e[5], 'scfend' => $e[7], 'ctgbeg' => $e[10], 'ctgend' => $e[12] );
				$HSP{$HSPcounter} = \%current_HSP;
			}
		}
		close ALIGNMENT;

		# inconsistent strand
		if($strand eq "inconsistent") {
			unlink("$outfolder/Alignments/alignments_$scaffold-$id.txt");
			next;
		}


		### Analyze HSP breaks
		for( my $i = 1; $i < $HSPcounter; $i++) {

			$total_aligned += $HSP{$i}{scfend} - $HSP{$i}{scfbeg} + 1;

			my $miss_scf = abs($HSP{$i+1}{scfbeg} - $HSP{$i}{scfend}) - 1;
			my $miss_ctg = abs($HSP{$i+1}{ctgbeg} - $HSP{$i}{ctgend}) - 1;

			if( $miss_scf >= 20 && $miss_ctg >= 20 ) {
				my $miss_scf_seq = substr($SEQ,  $HSP{$i}{scfend}, $miss_scf);
				my $N_scf_count  = $miss_scf_seq =~ s/([nN])/$1/gi;
				my $N_scf_perc   = $N_scf_count / length($miss_scf_seq);

				#my $miss_ctg_seq = substr($CTG{$id}, $HSP{$i}{ctgend}, $miss_ctg);
				#my $N_ctg_count  = $miss_ctg_seq =~ s/([nN])/$1/gi;
				#my $N_ctg_perc   = $N_ctg_count / length($miss_ctg_seq);

				my $ratio = min($miss_scf, $miss_ctg) / max($miss_scf, $miss_ctg);

				# Check if enough anchor sequence is alignable
				my $check_length = 1;
				if(
					(abs($HSP{$i}{scfend} - $HSP{$i}{scfbeg}) < 200) ||
					(abs($HSP{$i+1}{scfend} - $HSP{$i+1}{scfbeg}) < 200) ||
					(abs($HSP{$i}{ctgend} - $HSP{$i}{ctgbeg}) < 200) ||
					(abs($HSP{$i+1}{ctgend} - $HSP{$i+1}{ctgbeg}) < 200)
				){
					$check_length = 0;
				}
				
				# Check if a reasonable number of Ns can be corrected
				if( $N_scf_perc > 0.5 && $ratio > 0.5 && $check_length == 1) {

					my $plot_scf_seq = substr($SEQ,  $HSP{$i}{scfend} - 200, $miss_scf + 400);
					my $plot_ctg_seq = substr($CTG{$id}, $HSP{$i}{ctgend} - 200, $miss_ctg + 400);

					if($strand eq "-1") {
						print HDRFILE "Found it: $scaffold\t$id\n";
						$plot_ctg_seq = substr($CTG{$id}, $HSP{$i+1}{ctgbeg} - 200, $miss_ctg + 400);
						$plot_ctg_seq = revcompdna($plot_ctg_seq);
					}

					print HDRFILE "HDR\t$scaffold\t$id\t$i\t$strand\t" . $HSP{$i}{scfend} ."\t". $HSP{$i+1}{scfbeg} ."\t". 
						$HSP{$i}{ctgend} ."\t". $HSP{$i+1}{ctgbeg} . "\t$miss_scf\t$miss_ctg\t$plot_scf_seq\t$plot_ctg_seq\n";
				}
			}
		}


		if($total_aligned / length($CTG{$id}) > 0.95) {
			delete($CTG{$id});
		}
	}
}

close HDRFILE;

print STDERR "Finished Whole Genome Alignment\n";
exit(0);

sub revcompdna {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	#$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}
