#!/usr/bin/perl

use strict;
use warnings;


### User params
my $usage = "$0 ecotype assembly reference geneAnnotation chromosome maxDivBases outfolder\n";

my $ecotype   = shift or die $usage;
my $assembly  = shift or die $usage;
my $reference = shift or die $usage;
my $geneGFF   = shift or die $usage; 
my $chr       = shift or die $usage;
my $maxDivBp  = shift or die $usage;
my $outfolder = shift or die $usage;


### Init other variables
my $chr_int = substr($chr, 3);


### Create folders
mkdir "$outfolder" or die;
mkdir "$outfolder/Genes" or die;
mkdir "$outfolder/TEs" or die;
mkdir "$outfolder/Pseudogenes" or die;
mkdir "$outfolder/Alignments" or die;
mkdir "$outfolder/PartialGeneAlignments" or die;
mkdir "$outfolder/HSP_Alignments" or die;


### Run NUCMER
system("nucmer --mum -b $maxDivBp -g 90 -l 35 -c 80 -f --prefix=$outfolder/wgaMum $reference $assembly");

### Run delta filter
system("delta-filter -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");


### Parse reference file
open REFERENCE, $reference or die "Cannot open $reference\n";

my $id = "";
my $REF = "";

while(<REFERENCE>) {
        chomp;
	
	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		if($id eq $chr) {
			$REF .= $_;
		}
	}
}


### Parse assembly file
open ASSEMBLY, $assembly or die "Cannot open $assembly\n";

$id = "";
my %CTG = ();

while(<ASSEMBLY>) {
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
close ASSEMBLY;


### Parse gene annotation in GFF format
open GENE, $geneGFF or die "Cannot open $geneGFF\n";
my %geneAnn = ();
my %CDSAnn = ();
my %UTRAnn  = ();

while(<GENE>) {
	chomp;
	my @e = split("\t", $_);

	if( $e[0] eq $chr ) {

		# Gene
		if( $e[2] eq "gene" ) { 
			$geneAnn{$e[8]} = \@e;
		}

		# TE
		elsif( ($e[2] eq "transposable_element_gene") || ($e[2] eq "transposable_element") ) {
			$geneAnn{$e[8]} = \@e;
		}

		# Pseudogene
		if( $e[2] eq "pseudogene" ) {
			$geneAnn{$e[8]} = \@e;
		}

		# Exon
		elsif( $e[2] eq "CDS" ) {
			$CDSAnn{$e[8]} = \@e;
		}

		# UTR
		elsif( ($e[2] eq "five_prime_UTR") || ($e[2] eq "three_prime_UTR") ) {
			$UTRAnn{$e[8]} = \@e;
		}
	}
}


### Open result files and hashes
open SNPFILE, ">$outfolder/snp.txt" or die  "Cannot open $outfolder/snp.txt\n";
open DELFILE, ">$outfolder/del.txt" or die  "Cannot open $outfolder/del.txt\n";
open INSFILE, ">$outfolder/ins.txt" or die  "Cannot open $outfolder/ins.txt\n";
open HDRFILE, ">$outfolder/hdr.txt" or die  "Cannot open $outfolder/hdr.txt\n";
open AMBHDR, ">$outfolder/amb_hdr.txt" or die  "Cannot open $outfolder/amb_hdr.txt\n";
open AMBINS, ">$outfolder/amb_ins.txt" or die  "Cannot open $outfolder/amb_ins.txt\n";
open STUPIDO, ">$outfolder/stupido.txt" or die  "Cannot open $outfolder/stupido.txt\n";

my %SNP = ();
my %DEL = ();
my %INS = ();


### Create alignments for each scaffold ID and parse alignment results (e.g. SNPs, Indels, HDRs, gene annotation)
foreach $id (keys %CTG) {

	### Create and open alignment for one scaffold
	system("show-aligns -r $outfolder/wgaMum.filter $chr $id > $outfolder/Alignments/alignments_$id.txt");
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
			open FASTA, ">$outfolder/HSP_Alignments/$chr-$id-$HSPcounter.fa" or die "Cannot open $outfolder/HSP_Alignments/$chr-$id-$HSPcounter.fa\n";
			print FASTA ">". $chr ."-". $HSP{$HSPcounter}{refbeg} ."-". $HSP{$HSPcounter}{refend} ."\n$refseq\n";
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


			### Check if gene annotation is overlapping and print alignments of genes completely contained in one HSP
			foreach my $geneID (keys %geneAnn) {

				# Gene is completely covered by one HSP
				if( ($HSP{$HSPcounter}{refbeg} <= $geneAnn{$geneID}[3]) && ($geneAnn{$geneID}[4] <= $HSP{$HSPcounter}{refend}) ) {
					my @desc = split(";", $geneID);
					my ($junk, $name) = split("=", $desc[0]);

					# Walk to gene start and read gene alignment
					my $i = $HSP{$HSPcounter}{refbeg};
					my $gene_ref_aln = "";
					my $gene_ass_aln = "";
					my $anno_string  = "";

					for(my $j = 0; $j < length($refseq); $j++) {

						my $refnuc = substr($refseq, $j, 1);
						my $scaffnuc = substr($scaffseq, $j, 1);

						if( ($i >= $geneAnn{$geneID}[3]) && ($i <= $geneAnn{$geneID}[4]) ) {
							$gene_ref_aln .= $refnuc;
							$gene_ass_aln .= $scaffnuc;

							# Annotation string
							if(exists $CDSAnn{$i}) {
								$anno_string .= "E"
							}
							elsif(exists $UTRAnn{$i}) {
								$anno_string .= "U"
							}
							else {
								$anno_string .= "I"
							}
						}
						elsif($i > $geneAnn{$geneID}[4]) {
							last;
						}

						if($refnuc ne '.') { 
							$i++;
						}
					}

					if( $geneAnn{$geneID}[2] eq "gene" ) {
						open GENEALN, ">$outfolder/Genes/$name.fa" or die "Cannot open $outfolder/Genes/$name.fa";
						print GENEALN ">$name\n$gene_ref_aln\n>$id\n$gene_ass_aln\n>TAIR8\n$anno_string\n";
						close GENEALN;
					}
					elsif( ($geneAnn{$geneID}[2] eq "transposable_element_gene") || ($geneAnn{$geneID}[2] eq "transposable_element") ) {
						open GENEALN, ">$outfolder/TEs/$name.fa" or die "Cannot open $outfolder/TEs/$name.fa";
						print GENEALN ">$name\n$gene_ref_aln\n>$id\n$gene_ass_aln\n>TAIR8\n$anno_string\n";
						close GENEALN;
					}
					elsif( $geneAnn{$geneID}[2] eq "pseudogene" ) {
						open GENEALN, ">$outfolder/Pseudogenes/$name.fa" or die "Cannot open $outfolder/Pseudogenes/$name.fa";
						print GENEALN ">$name\n$gene_ref_aln\n>$id\n$gene_ass_aln\n>TAIR8\n$anno_string\n";
						close GENEALN;
					}
				}

				# Gene overlaps with one or more HSPs
				elsif (	
					( ($geneAnn{$geneID}[3] <= $HSP{$HSPcounter}{refend}) && ($HSP{$HSPcounter}{refend} <= $geneAnn{$geneID}[4]) ) ||
					( ($geneAnn{$geneID}[3] <= $HSP{$HSPcounter}{refbeg}) && ($HSP{$HSPcounter}{refbeg} <= $geneAnn{$geneID}[4]) )
				) {
					my @desc = split(";", $geneID);
					my ($junk, $name) = split("=", $desc[0]);

					# Walk to gene start and read gene alignment
					my $i = $HSP{$HSPcounter}{refbeg};
					my $gene_ref_aln = "";
					my $gene_ass_aln = "";

					for(my $j = 0; $j < length($refseq); $j++) {
						
						my $refnuc = substr($refseq, $j, 1);
						my $scaffnuc = substr($scaffseq, $j, 1);

						if( ($i >= $geneAnn{$geneID}[3]) && ($i <= $geneAnn{$geneID}[4]) ) {
							$gene_ref_aln .= $refnuc;
							$gene_ass_aln .= $scaffnuc;
						}
						elsif($i > $geneAnn{$geneID}[4]) {
							last;
						}
						
						if($refnuc ne '.') {
							$i++;
						}
					}

					if(-e "$outfolder/PartialGeneAlignments/$name.fa" ) {
						
						# Print overlapping scaffold sequence
						open GENEALN, ">>$outfolder/PartialGeneAlignments/$name.fa" or die "Cannot open $outfolder/PartialGeneAlignments/$name.fa";
						print GENEALN ">$id HSP-$HSPcounter\n$gene_ass_aln\n";
						close GENEALN;
					}
					else {
						# Print full gene sequence and overlapping scaffold sequence
						my $gene_seq = substr( $REF, ($geneAnn{$geneID}[3] - 1), ($geneAnn{$geneID}[4] - $geneAnn{$geneID}[3] + 1) );  
						open GENEALN, ">$outfolder/PartialGeneAlignments/$name.fa" or die "Cannot open $outfolder/PartialGeneAlignments/$name.fa";
						print GENEALN ">$name\n$gene_seq\n>$id HSP-$HSPcounter\n$gene_ass_aln\n";
						close GENEALN;
					}
				}
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
		if( ($miss_ref <= 100000) && ($miss_scaff <= 100000) ) {

			# Reference and Scaffold have additional sequence missing between the HSPs
			if(  ($miss_ref >= 20) && ($miss_scaff >= 20) ) {

				my $miss_ref_seq  = substr($REF,  $HSP{$i}{refend}, $miss_ref);
				my $miss_scaff_seq = substr($CTG{$id}, $HSP{$i}{scaffend}, $miss_scaff);
				
				my $N_ref_count = $miss_ref_seq =~ s/([nN])/$1/gi;
				my $N_ref_perc  = $N_ref_count / length($miss_ref_seq);

				my $N_scaff_count = $miss_scaff_seq =~ s/([nN])/$1/gi;
				my $N_scaff_perc  = $N_scaff_count / length($miss_scaff_seq);

				if( ($N_ref_perc <= 0.05) && ($N_scaff_perc<= 0.05) ) {
					print HDRFILE "$ecotype\tHDR\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t". 
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t$miss_scaff_seq\n";
				}
				else {
					print AMBHDR "$ecotype\tHDR\t$chr_int\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t$id\t". 
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t$miss_scaff_seq\n";
				}
			}
	
			# Misalignments with negative HSP distance
			elsif( ($miss_ref < -20) || ($miss_scaff < -20) ) {
				print STUPIDO "$ecotype\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
			}

			# Reference has additional sequence missing between the HSPs (equals a deletion)
			elsif($miss_ref >= 20) {

				my $miss_ref_seq  = substr($REF,  $HSP{$i}{refend}, $miss_ref);

				my $N_count = $miss_ref_seq =~ s/([nN])/$1/gi;
				my $N_perc  = $N_count / length($miss_ref_seq);

				if( $N_perc <= 0.05) {
					print HDRFILE "$ecotype\tDEL\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t#\n";
				}
				else {
					print AMBHDR "$ecotype\tDEL\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t$miss_ref_seq\t#\n";
				}
			}

			# Scaffold has additional sequence missing between the HSPs (equals an insertion)
			elsif($miss_scaff >= 20) {

				my $miss_scaff_seq = substr($CTG{$id}, $HSP{$i}{scaffend}, $miss_scaff);

				my $N_count = $miss_scaff_seq =~ s/([nN])/$1/gi;
				my $N_perc  = $N_count / length($miss_scaff_seq);
				
				if( $N_perc <= 0.05) {
					print HDRFILE "$ecotype\tINS\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t#\t$miss_scaff_seq\n";
				}
				else {
					print AMBHDR "$ecotype\tINS\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
						$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\t#\t$miss_scaff_seq\n";
				}
			}

			# Scaffold and Reference miss less then 20bp - why HSP break?
			else {
				print STUPIDO "$ecotype\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
					$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
			}
		}

		# Huge HSP breaks
		else {
			print STUPIDO "$ecotype\t$chr_int\t$id\t$i\t" . $HSP{$i}{refend} ."\t". $HSP{$i+1}{refbeg} ."\t".
				$HSP{$i}{scaffend} ."\t". $HSP{$i+1}{scaffbeg} . "\t$miss_ref\t$miss_scaff\n";
		}
	}

	close ALIGNMENT;
}



### Sort, concatenate and print SNPs and Indels

# SNPs
foreach my $pos ( sort {$a<=>$b} keys %SNP) {
	print SNPFILE "$ecotype\t$chr_int\t$pos\t" . $SNP{$pos}{id} ."\t". $SNP{$pos}{position} ."\t". $SNP{$pos}{nuc} . "\n";
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
			print DELFILE "$ecotype\t$chr_int\t$delbeg\t$delend\t" . $DEL{$pos}{id} ."\t". $DEL{$pos}{position} . "\t$scaff_end\t$len\t$delseq\n";
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

	#if($INS{$pos}{seq} !~ /[nN]/) {
	if( $N_perc <= 0.05) {
		print INSFILE "$ecotype\t$chr_int\t$pos\t$end\t" . $INS{$pos}{id} . "\t$scaff_begin\t" . $INS{$pos}{position} . "\t$len\t" . $INS{$pos}{seq} . "\n";
	}
	else {
		print AMBINS "$ecotype\tINS\t$chr_int\t$pos\t$end\t" . $INS{$pos}{id} . "\t$scaff_begin\t" . $INS{$pos}{position} . "\t$len\t" . $INS{$pos}{seq} . "\n";
	}
}

close SNPFILE;
close DELFILE;
close INSFILE;
close HDRFILE;
close AMBHDR;
close AMBINS;

print STDERR "Finished Whole Genome Alignment\n";
exit(0);
