#!/usr/bin/perl -w

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

use SNPlist;
use IndelList;
use HDRList;
use GeneSNPlist;


### Variables
my $refseq_file = "";
my $gff         = "";
my $snp_file    = "";
my $del_file    = "";
my $ins_file    = "";
my $hdr_file    = "";
my $chromosome;



### Get command line options
my %CMD;
GetCom();

print STDERR "Reading input files\n";

### Get SNP lists
my $snps = new SNPlist();
$snps->get($chromosome, $snp_file);

### Get deletion list
my $deletions = new IndelList();
if($del_file ne "") { $deletions->get($chromosome, $del_file); }

### Get insertion list
my $insertions = new IndelList();
if($ins_file ne "") { $insertions->get($chromosome, $ins_file); }

### Get HDR list
my $HDRs = new HDRList();
if($hdr_file ne "") { $HDRs->get($chromosome, $hdr_file); }


### Functional analysis of SNPs and indels
my %coding_ann = ();
my %gene_ann = ();
my %seq_type = ();

if( ($gff ne "") && ($refseq_file ne "") ) {

	### Get gene annotation from gff
	open GFF, $gff or die "Cannot open gff file\n";
	while( <GFF> ) {

		# GFF format: chr, source, seq_type, start, end, score, orientation, frame, description
		my @columns = split(" ", $_);

		# Quick fix for A.thaliana TAIR8 annotation, disable for other species
		if($columns[0] =~ /Chr/) { $columns[0] = substr($columns[0], 3); }

		# Check if gene is on selected chromosome
		if( $columns[0] eq $chromosome ) {

			# Split decription column
			my @desc = split(";", $columns[8]);

			# Get gene locus
			if( $columns[2] eq "gene" ) {
				if( substr($desc[0], 0, 2) eq "ID") {
					my ($junk, $gene_name)  = split("=", $desc[0]);
					if(! exists $coding_ann{$gene_name}) {
						#my %gene = ();
						#$coding_ann{$gene_name} = \%gene;

						my @gene_locus = ($columns[3], $columns[4], $columns[6], "");
						$gene_ann{$gene_name} = \@gene_locus;

						my %type = ();
						$seq_type{$gene_name} = \%type;
						$seq_type{$gene_name}{"gene"} = $columns[3] . "-" . $columns[4];
					}
				}
			}

			# Get mRNA locus
			elsif( $columns[2] eq "mRNA" ) {
				# TODO, currently not needed
			}

			# Get coding sequence
			elsif( $columns[2] eq "CDS" ) {
				my ($parent, $protein_name) = split(",", $desc[0]);
				my ($junk, $gene_id) = split("=", $parent);
				my $isoform = 1;
				my $gene_name = $gene_id;

				if( $gene_id =~ /\./ ) {
					($gene_name, $isoform) = split(/\./, $gene_id);
				}

				if($isoform == 1) {
					if(! exists $coding_ann{$gene_name}) {
						my %gene = ();
						$coding_ann{$gene_name} = \%gene;
					}

					my @cds = ($columns[3], $columns[4], $columns[6], "");
					$coding_ann{$gene_name}{$columns[3]} = \@cds;

					$seq_type{$gene_name}{"CDS"} = $columns[3] . "-" . $columns[4];
				}
			}

			# Get other sequence types
			elsif( ($columns[2] eq "exon") || ($columns[2] eq "five_prime_UTR") || ($columns[2] eq "three_prime_UTR") ) {
				my ($junk, $gene_id) = split("=", $desc[0]);
				my $isoform = 1;
				my $gene_name = $gene_id;
				if( $gene_id =~ /\./ ) {
					($gene_name, $isoform) = split(/\./, $gene_id);
				}

				if($isoform == 1) {
					$seq_type{$gene_name}{$columns[2]} = $columns[3] . "-" . $columns[4];
				}
			}
		}
	}


	### Load chromosome
	my $chr_seq = "";
	open GENOME, $refseq_file or die "Cannot open reference sequence file\n";
	while( <GENOME> ) {
		chomp;
		if(substr($_, 0, 1)  eq ">") {
			my $current_chr = substr($_, 1);
			$current_chr =~ s/\s+//g;

			if($current_chr eq $chromosome) {
				while(<GENOME>) {
					chomp;
					if(substr($_, 0, 1)  eq ">") { 
						last; 
					}
					else {
						$chr_seq .= $_;
					}
				}
			}
		}

		if($chr_seq ne "" ) { last; }
	}


	### Extract gene seq
	foreach my $gene_name ( sort keys %coding_ann ) {
		$gene_ann{$gene_name}[3] = substr($chr_seq, $gene_ann{$gene_name}[0] - 1, $gene_ann{$gene_name}[1] - $gene_ann{$gene_name}[0] + 1);

		foreach my $start ( sort {$a <=> $b} keys %{$coding_ann{$gene_name}} ) {
			$coding_ann{$gene_name}{$start}[3] = substr($chr_seq, $start - 1, $coding_ann{$gene_name}{$start}[1] - $start + 1);
		}
	}

	### Add gene SNP and Indel annotation
	foreach my $gene_name ( sort keys %gene_ann ) {
		
		# Coding SNP protein change
		if( exists $coding_ann{$gene_name} ) {
			my $gene = new GeneSNPlist($snps);
			my $results = $gene->get_gene_snps( $chromosome, $gene_ann{$gene_name}[0], $gene_ann{$gene_name}[1],
							$gene_ann{$gene_name}[2], $gene_name, 1, $gene_ann{$gene_name}[3], %{$coding_ann{$gene_name}});
			$gene->get_protein_changes();

		}

		# Calculate SNP and Indel annotation type
		foreach my $seq_type ( keys %{$seq_type{$gene_name}} ) {
			my ($seq_type_start, $seq_type_end) = split( "-", $seq_type{$gene_name}{$seq_type} );

			# SNPs
			foreach my $snp_pos ( sort {$a <=> $b} keys %{$snps->{snps}} ) {

				if( ($snp_pos >= $seq_type_start) && ($snp_pos <= $seq_type_end) ) {

					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$snps->{snps}{$snp_pos}{stype} = $seq_type;
					}

					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $snps->{snps}{$snp_pos}{stype} ne "CDS") {
							$snps->{snps}{$snp_pos}{stype} = $seq_type;
							$snps->{snps}{$snp_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic
					elsif( $seq_type eq "gene" ) {
						if(	($snps->{snps}{$snp_pos}{stype} ne "CDS") && 
							($snps->{snps}{$snp_pos}{stype} ne "five_prime_UTR") && 
							($snps->{snps}{$snp_pos}{stype} ne "three_prime_UTR")
						){
							$snps->{snps}{$snp_pos}{stype} = "intronic";
							$snps->{snps}{$snp_pos}{gene_id} = $gene_name;
						}
					}
				}
			}

			# Deletions
			foreach my $del_pos ( sort {$a <=> $b} keys %{$deletions->{indels}} ) {
				if( 
					( ($del_pos >= $seq_type_start) && ($del_pos <= $seq_type_end) ) ||
					( ($deletions->{indels}{$del_pos}{end} >= $seq_type_start) && ($deletions->{indels}{$del_pos}{end} <= $seq_type_end) )
				) {
					
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$deletions->{indels}{$del_pos}{stype} = $seq_type;
						$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
					}

					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $deletions->{indels}{$del_pos}{stype} ne "CDS") {
							$deletions->{indels}{$del_pos}{stype} = $seq_type;
							$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic
					elsif( $seq_type eq "gene" ) {
						if(     ($deletions->{indels}{$del_pos}{stype} ne "CDS") &&
							($deletions->{indels}{$del_pos}{stype} ne "five_prime_UTR") &&
							($deletions->{indels}{$del_pos}{stype} ne "three_prime_UTR")
						){
							$deletions->{indels}{$del_pos}{stype} = "intronic";
							$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
						}
					}
				}
			}

			# Insertions
			foreach my $ins_pos ( sort {$a <=> $b} keys %{$insertions->{indels}} ) {
				
				if( 
					( ($ins_pos >= $seq_type_start) && ($ins_pos <= $seq_type_end) ) ||
					( ($insertions->{indels}{$ins_pos}{end} >= $seq_type_start) && ($insertions->{indels}{$ins_pos}{end} <= $seq_type_end) )
				) {
					
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$insertions->{indels}{$ins_pos}{stype} = $seq_type;
						$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
					}


					# Set sequence type to UTR
					elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
						if( $insertions->{indels}{$ins_pos}{stype} ne "CDS") {
							$insertions->{indels}{$ins_pos}{stype} = $seq_type;
							$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
						}
					}

					# Set sequence type to intronic
					elsif( $seq_type eq "gene" ) {
						if(	($insertions->{indels}{$ins_pos}{stype} ne "CDS") &&
							($insertions->{indels}{$ins_pos}{stype} ne "five_prime_UTR") &&
							($insertions->{indels}{$ins_pos}{stype} ne "three_prime_UTR")
						){
							$insertions->{indels}{$ins_pos}{stype} = "intronic";
							$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
						}
					}
				}
			}

			# HDR
			foreach my $hdr_beg ( sort {$a <=> $b} keys %{$HDRs->{indels}} ) {
				my $hdr_end = $HDRs->{indels}{$hdr_beg}{end};

				for(my $i = $hdr_beg; $i <= $hdr_end; $i++) {
					
					if( ($i >= $seq_type_start) && ($i <= $seq_type_end) ) {
						
						# Set sequence type to CDS
						if( $seq_type eq "CDS" ) {
							$HDRs->{indels}{$hdr_beg}{stype} = "CDS";
							$HDRs->{indels}{$hdr_beg}{cds_absence}{$gene_name} = 1;
						}

						# Set sequence type to UTR
						elsif( ($seq_type eq "five_prime_UTR") || ($seq_type eq "three_prime_UTR") ) {
							if( $HDRs->{indels}{$hdr_beg}{stype} ne "CDS") {
								$HDRs->{indels}{$hdr_beg}{stype} = $seq_type;
							}
							$HDRs->{indels}{$hdr_beg}{utr_absence}{$gene_name} = 1;
						}

						# Set sequence type to intronic
						elsif( $seq_type eq "gene" ) {
							if(	($HDRs->{indels}{$hdr_beg}{stype} ne "CDS") &&
								($HDRs->{indels}{$hdr_beg}{stype} ne "five_prime_UTR") &&
								($HDRs->{indels}{$hdr_beg}{stype} ne "three_prime_UTR")
							){
								$HDRs->{indels}{$hdr_beg}{stype} = "intronic";
							}
							$HDRs->{indels}{$hdr_beg}{gene_absence}{$gene_name} = 1;
						}
					}
				}
			}
		}
	}
}


### Print SNP results
open SNPOUT, ">snp.annotation.txt" or die "Cannot open SNP output file\n";

foreach my $pos (sort {$snps->{snps}{$a}{position} <=> $snps->{snps}{$b}{position}} keys %{$snps->{snps}} ) {
	print SNPOUT "$chromosome\t$pos\t" . $snps->{snps}{$pos}{ref_base} ."\t". $snps->{snps}{$pos}{new_base}; 

	if( ($gff ne "") && ($refseq_file ne "") ) {

		print SNPOUT "\t" . $snps->{snps}{$pos}{stype};

		if($snps->{snps}{$pos}{gene_id} ne "NA") {
			print SNPOUT "\t" . $snps->{snps}{$pos}{gene_id} . "\t1";
		}

		if($snps->{snps}{$pos}{cds_pos} != 0) {
			my $syn_nonsyn = "Syn";
			if($snps->{snps}{$pos}{ns_change} == 1) { $syn_nonsyn = "Nonsyn"; }

			print SNPOUT "\t" . $snps->{snps}{$pos}{cds_pos} . "\t" . $snps->{snps}{$pos}{codon_pos} . "\t$syn_nonsyn\t" . 
				$snps->{snps}{$pos}{ref_aa} . "\t" . $snps->{snps}{$pos}{new_aa};
		}
	}

	print SNPOUT "\n";
}
close SNPOUT;


### Print deletion results
if($del_file ne "") {
	open DELOUT, ">del.annotation.txt" or die "Cannot open deletion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$deletions->{indels}} ) {
		
		print DELOUT "$chromosome\t$start\t" . $deletions->{indels}{$start}{end} ."\t". $deletions->{indels}{$start}{seq}; 
				
		if( ($gff ne "") && ($refseq_file ne "") ) {	
			print DELOUT "\t". $deletions->{indels}{$start}{stype};

			if( $deletions->{indels}{$start}{gene_id} ne "NA") {
				print DELOUT "\t". $deletions->{indels}{$start}{gene_id};
			}
		}
	
		print DELOUT "\n";
	}

	close DELOUT;
}


### Print insertion results
if($ins_file ne "") {
	open INSOUT, ">ins.annotation.txt" or die "Cannot open insertion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$insertions->{indels}} ) {

		print INSOUT "$chromosome\t$start\t" . $insertions->{indels}{$start}{end} ."\t". $insertions->{indels}{$start}{seq};
				
		if( ($gff ne "") && ($refseq_file ne "") ) {
			print INSOUT "\t". $insertions->{indels}{$start}{stype};

			if( $insertions->{indels}{$start}{gene_id} ne "NA") {
				print INSOUT "\t". $insertions->{indels}{$start}{gene_id};
			}
		}

		print INSOUT "\n";
	}
	
	close INSOUT;
}

### Print HDR results
if($hdr_file ne "") {
	open HDROUT, ">hdr.annotation.txt" or die "Cannot open HDR output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$HDRs->{indels}} ) {

		print HDROUT "$chromosome\t$start\t" . $HDRs->{indels}{$start}{end} . "\t" . $HDRs->{indels}{$start}{dtype};

		if( ($gff ne "") && ($refseq_file ne "") ) {
		
			print HDROUT "\t". $HDRs->{indels}{$start}{stype};

			my $gene_list = "";
			my $utr_list = "";
			my $cds_list = "";

			foreach my $gene_id ( sort keys %{$HDRs->{indels}{$start}{gene_absence}} ) {
				$gene_list .= $gene_id . ",";
			}
			foreach my $gene_id ( sort keys %{$HDRs->{indels}{$start}{utr_absence}} ) {
				$utr_list .= $gene_id . ",";
			}
			foreach my $gene_id ( sort keys %{$HDRs->{indels}{$start}{cds_absence}} ) {
				$cds_list .= $gene_id . ",";
			}

			if($gene_list eq "") { $gene_list = "NA"; }
			else { chop($gene_list); }

			if($utr_list eq "") { $utr_list = "NA"; }
			else { chop($utr_list); }

			if($cds_list eq "") { $cds_list = "NA"; }
			else { chop($cds_list); }

			print HDROUT "\t$gene_list\t$utr_list\t$cds_list";
		}

		print HDROUT "\n";
	}
	close HDROUT;
}

exit(0);


### Read command line parameters
sub GetCom {
  my @usage = ("\nUsage: $0

Variation files (mandatory):
--snp      STRING     SNP file in WGA format (see Mummer_Parser.pl)
--del      STRING     Deletion file
--ins      STRING     Insertion file
--hdr      STRING     HDR file

Functional SNP and indel annotation (mandatory):
--chrom    STRING     Chromosome
--genome   STRING     Reference sequence file (chromosome names have to be equal to SNP file)
--gff      STRING     Gene annotation in GFF format

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "snp=s", "del=s", "ins=s", "hdr=s", "chrom=s", "genome=s", "gff=s");

	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify del file\n") unless defined($CMD{del});
	die("Please specify ins file\n") unless defined($CMD{ins});
	die("Please specify hdr file\n") unless defined($CMD{hdr});
	die("Please specify chromosome of target region\n") unless defined($CMD{chrom});

	$snp_file    = $CMD{snp};
	$del_file    = $CMD{del};
	$ins_file    = $CMD{ins};
	$hdr_file    = $CMD{hdr};
	$chromosome  = $CMD{chrom};
	$refseq_file = $CMD{genome};
	$gff         = $CMD{gff};
}

