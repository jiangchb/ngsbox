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




### User params
my $usage = "\n$0 contig_file scaffold_file scaffold maxDivBases outfolder\n\n";

my $contigs   = shift or die $usage;
my $scaffolds = shift or die $usage;
my $scaffold  = shift or die $usage;
my $maxDivBp  = shift or die $usage;
my $outfolder = shift or die $usage;



### Create folders
mkdir "$outfolder" or die;
mkdir "$outfolder/Alignments" or die;
mkdir "$outfolder/HSP_Alignments" or die;


### Run NUCMER
system("nucmer --mum -b $maxDivBp -g 90 -l 35 -c 80 -f --prefix=$outfolder/wgaMum $scaffolds $contigs");

### Run delta filter
system("delta-filter -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");


### Parse scaffolds file
open SCAFFOLDS, $scaffolds or die "Cannot open $scaffolds\n";

my $id = "";
my $SEQ = "";
while(<SCAFFOLDS>) {
        chomp;
	
	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		if($id eq $scaffold) {
			$SEQ .= $_;
		}
	}
}


### Parse contigs file
open CONTIGS, $contigs or die "Cannot open $contigs\n";

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
open SNPFILE, ">$outfolder/snp.txt" or die  "Cannot open $outfolder/snp.txt\n";
open DELFILE, ">$outfolder/del.txt" or die  "Cannot open $outfolder/del.txt\n";
open INSFILE, ">$outfolder/ins.txt" or die  "Cannot open $outfolder/ins.txt\n";
open HDRFILE, ">$outfolder/hdr.txt" or die  "Cannot open $outfolder/hdr.txt\n";
open STUPIDO, ">$outfolder/stupido.txt" or die  "Cannot open $outfolder/stupido.txt\n";

my %SNP = ();
my %DEL = ();
my %INS = ();


### Create alignments for each scaffold ID and parse alignment results (e.g. SNPs, Indels, HDRs, gene annotation)
foreach $id (keys %CTG) {

	# TEST
	#system("show-aligns -r $outfolder/wgaMum.filter $scaffold $id > $outfolder/Alignments/alignments_$id.txt") == 0 or print "system failed: $?\n";

	### Create and open alignment for one scaffold
	my $exit_code = system("show-aligns -r $outfolder/wgaMum.filter $scaffold $id > $outfolder/Alignments/alignments_$id.txt");
	if($exit_code != 0) {
		unlink("$outfolder/Alignments/alignments_$id.txt");
		next;
	}

	open ALIGNMENT, "$outfolder/Alignments/alignments_$id.txt" or die "Cannot open $outfolder/Alignments/alignments_$id.txt\n";


	### Init variables
	my $state = 0;
	my $junk_line = "";
	my $refseq = "";
	my $scaffseq = "";
	my $HSPcounter = 0;
	my %HSP = ();


	### Parse alignment file of one scaffold
	while(<ALIGNMENT>) {
		chomp;

		### New HSP found ----------------------------------------------------------------------
		if( $_ =~ /BEGIN/ ) {
			$HSPcounter++;
			my @e = split(" ", $_);
			my %current_HSP = ( 'refbeg' => $e[5], 'refend' => $e[7], 'scaffbeg' => $e[10], 'scaffend' => $e[12] );
			$HSP{$HSPcounter} = \%current_HSP;

			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1;
		}


		### HSP finished, process information --------------------------------------------------
		elsif(  $_ =~ /END/ ) {

			### Print Alignment in Multi-Fasta format
			open FASTA, ">$outfolder/HSP_Alignments/$scaffold-$id-$HSPcounter.fa" or die "Cannot open $outfolder/HSP_Alignments/$scaffold-$id-$HSPcounter.fa\n";
			print FASTA ">". $scaffold ."-". $HSP{$HSPcounter}{refbeg} ."-". $HSP{$HSPcounter}{refend} ."\n$refseq\n";
			print FASTA ">". $id ."-". $HSP{$HSPcounter}{scaffbeg} ."-". $HSP{$HSPcounter}{scaffend} ."\n$scaffseq\n"; 
			close FASTA;

			
			### Variants (SNPs, Indels)
			my $r = $HSP{$HSPcounter}{refbeg};
			my $s = $HSP{$HSPcounter}{scaffbeg};
			for( my $j = 0; $j < length($refseq); $j++ ) {

				my $refnuc = substr($refseq, $j, 1);
				my $scaffnuc = substr($scaffseq, $j, 1);

				if( $refnuc ne $scaffnuc ) {

					if( ($refnuc ne '.') && ($scaffnuc ne '.') && ($scaffnuc ne 'n')) {
						$SNP{$r}{id} = $id;
						$SNP{$r}{position} = $s;
						$SNP{$r}{nuc} = $refnuc ."\t". $scaffnuc;
					}
					elsif( $refnuc eq '.') {
						$INS{$r}{id} = $id;
						$INS{$r}{position} = $s;
						$INS{$r}{seq} .= $scaffnuc;
					}
					elsif( $scaffnuc eq '.') {
						$DEL{$r}{id} = $id;
						$DEL{$r}{position} = $s;
						$DEL{$r}{nuc} = $refnuc;
					}
				}

				if($refnuc ne '.') { $r++; }
				if($scaffnuc ne '.') { $s++; }
			}


			### Reset
			$refseq = "";
			$scaffseq = "";
			$state = 0;
		}


		### Correct for stupid extra empty line after HSP-END ----------------------------------
		elsif( $_ eq "" ) {
			# Do nothing
		}


		### Reference line of alignment --------------------------------------------------------
		elsif( $state == 1 ) {
			my ($refpos, $seq) = split(" ", $_);
			$refseq .= $seq;
			$state = 2;
		}


		###Scaffold line of alignment ----------------------------------------------------------
		elsif( $state == 2 ) { 
			my ($scaffpos, $seq) = split(" ", $_);
			$scaffseq .= $seq;
			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1; 
		}
	}


	### Analyze HSP breaks
	for( my $i = 1; $i < $HSPcounter; $i++) {

		my $miss_ref = $HSP{$i+1}{refbeg} - $HSP{$i}{refend} - 1;
		my $miss_scaff = $HSP{$i+1}{scaffbeg} - $HSP{$i}{scaffend} - 1; 

		# HSP break of acceptable size
		if( ($miss_ref <= 30000) && ($miss_scaff <= 30000) ) {

			# Reference and Scaffold have additional sequence missing between the HSPs
			if(  ($miss_ref >= 10) && ($miss_scaff >= 10) ) {

				my $miss_ref_seq  = substr($SEQ,  $HSP{$i}{refend}, $miss_ref);
				my $miss_scaff_seq = substr($CTG{$id}, $HSP{$i}{scaffend}, $miss_scaff);
				
				my $N_ref_count = $miss_ref_seq =~ s/([nN])/$1/gi;
				my $N_ref_perc  = $N_ref_count / length($miss_ref_seq);

				my $N_scaff_count = $miss_scaff_seq =~ s/([nN])/$1/gi;
				my $N_scaff_perc  = $N_scaff_count / length($miss_scaff_seq);

				print HDRFILE "HDR\t$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t". 
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t$miss_scaff_seq\n";
			}
	
			# Misalignments with negative HSP distance
			elsif( ($miss_ref < -20) || ($miss_scaff < -20) ) {
				print STUPIDO "$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
			}

			# Reference has additional sequence missing between the HSPs (equals a deletion)
			elsif($miss_ref >= 10) {

				my $miss_ref_seq  = substr($SEQ,  $HSP{$i}{refend}, $miss_ref);

				my $N_count = $miss_ref_seq =~ s/([nN])/$1/gi;
				my $N_perc  = $N_count / length($miss_ref_seq);

				print HDRFILE "DEL\t$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t#\n";
			}

			# Scaffold has additional sequence missing between the HSPs (equals an insertion)
			elsif($miss_scaff >= 10) {

				my $miss_scaff_seq = substr($CTG{$id}, $HSP{$i}{scaffend}, $miss_scaff);

				my $N_count = $miss_scaff_seq =~ s/([nN])/$1/gi;
				my $N_perc  = $N_count / length($miss_scaff_seq);
				
				print HDRFILE "INS\t$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t#\t$miss_scaff_seq\n";
			}

			# Scaffold and Reference miss less then 10bp - why HSP break?
			else {
				print STUPIDO "$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
			}
		}

		# Huge HSP breaks
		else {
			print STUPIDO "$scaffold\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
				$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
		}
	}

	close ALIGNMENT;
}



### Sort, concatenate and print SNPs and Indels

# SNPs
foreach my $pos ( sort {$a<=>$b} keys %SNP) {
	print SNPFILE "$scaffold\t$pos\t" . $SNP{$pos}{id} ."\t". $SNP{$pos}{position} ."\t". $SNP{$pos}{nuc} . "\n";
}

# Deletions
my $delbeg = -1;
my $delend = -1;
my $delseq = "";

foreach my $pos ( sort {$a<=>$b} keys %DEL) {
	if( $pos != ($delend + 1) ) {
		if($delend != -1) {
			my $len = length($delseq);
			my $scaff_end = $DEL{$pos}{position} + 1;
			print DELFILE "$scaffold\t$delbeg\t$delend\t" . $DEL{$pos}{id} ."\t". $DEL{$pos}{position} . "\t$scaff_end\t$len\t$delseq\n";
		}
		$delbeg = $pos;
		$delseq = "";
	}
	$delend = $pos;
	$delseq .= $DEL{$pos}{nuc};
}

# Insertions
foreach my $pos ( sort {$a<=>$b} keys %INS) {
	my $end = $pos + 1;
	my $len = length($INS{$pos}{seq});
	my $scaff_begin = $INS{$pos}{position} - $len + 1;

	my $N_count = $INS{$pos}{seq} =~ s/([nN])/$1/gi;
	my $N_perc  = $N_count / length($INS{$pos}{seq});

	print INSFILE "$scaffold\t$pos\t$end\t" . $INS{$pos}{id} . "\t$scaff_begin\t" . $INS{$pos}{position} . "\t$len\t" . $INS{$pos}{seq} . "\n";
}

close SNPFILE;
close DELFILE;
close INSFILE;
close HDRFILE;

print STDERR "Finished Whole Genome Alignment\n";
exit(0);
