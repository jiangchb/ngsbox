#!/usr/bin/perl -w

# --------------------------------------------------------------------
# NGSBox: Variant Annotation
#
# Annotate SNPs, indels, SVs and CNVs against any eukaryotic gene
# annotation in GFF format
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

### TODO adapt for different species, most importantly for Ath

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
my $format      = "";
my $snp_file    = "";
my $del_file    = "";
my $ins_file    = "";
my $hdr_file    = "";
my $outfolder   = "";
my $chromosome;

print STDERR "Reading command line parameters\n";

### Get command line options
my %CMD;
GetCom();

### Get SNP lists
print STDERR "Reading SNP file\n";
my $snps = new SNPlist();
$snps->get($chromosome, $snp_file, $format);


### Get deletion list
print STDERR "Reading deletion file\n";
my $deletions = new IndelList();
if($del_file ne "") { $deletions->get($chromosome, $del_file, $format); }


### Get insertion list
print STDERR "Reading insertion file\n";
my $insertions = new IndelList();
if($ins_file ne "") { $insertions->get($chromosome, $ins_file, $format); }


### Get HDR list
print STDERR "Reading HDR file\n";
my $HDRs = new HDRList();
if($hdr_file ne "") { $HDRs->get($chromosome, $hdr_file, $format); }


### Functional analysis of SNPs and indels
my %coding_ann = ();
my %gene_ann = ();
my %seq_type = ();


print STDERR "Reading GFF file\n";

### Get gene annotation from gff
open GFF, $gff or die "Cannot open gff file\n";
while( <GFF> ) {

	# GFF format: chr, source, seq_type, start, end, score, orientation, frame, description
	my @columns = split("\t", $_);

	# Translate to shore chromosome format: remove "chr"
	if($columns[0] =~ /_/) { next; }
	elsif($columns[0] =~ /chr\d/) { $columns[0] = substr($columns[0], 3); }
	elsif($columns[0] =~ /chrX/)  { $columns[0] = 23; }
	elsif($columns[0] =~ /chrY/)  { $columns[0] = 24; }
	elsif($columns[0] =~ /chrM/)  { $columns[0] = 25; }


	# Check if gene is on selected chromosome
	if( $columns[0] eq $chromosome ) {

		# Split decription column
		my @desc = split(";", $columns[8]);
		my ($junk, $gene_name)  = split(" ", $desc[0]);
		$gene_name =~ s/"//g;

		# Add new gene
		if(! exists $gene_ann{$gene_name}) {
			my @gene_locus = ($columns[3], $columns[4], $columns[6], "");
			$gene_ann{$gene_name} = \@gene_locus;

			my %type = ();
			$seq_type{$gene_name} = \%type;
			my @seq_type_locus = ($columns[3], $columns[4]);
			$seq_type{$gene_name}{"gene"} = \@seq_type_locus;
		}


		# Start codon
		if( $columns[2] eq "start_codon" ) {
			# NIX
		}

		# Stop codon
		elsif( $columns[2] eq "stop_codon" ) {
			# NIX
		}

		# Get coding sequence
		elsif( $columns[2] eq "CDS" ) {
			my $isoform = 1;

			if( $gene_name =~ /\./ ) {
				($gene_name, $isoform) = split(/\./, $gene_name);
			}

			if($isoform == 1) {
				if(! exists $coding_ann{$gene_name}) {
					my %gene = ();
					$coding_ann{$gene_name} = \%gene;
				}

				my @cds = ($columns[3], $columns[4], $columns[6], "");
				$coding_ann{$gene_name}{$columns[3]} = \@cds;

				my @seq_type_locus = ($columns[3], $columns[4]);
				$seq_type{$gene_name}{"CDS"} = \@seq_type_locus;
			}
		}

		# Get exons sequence types
		elsif( $columns[2] eq "exon" ) {
			my $isoform = 1;
			if( $gene_name =~ /\./ ) {
				($gene_name, $isoform) = split(/\./, $gene_name);
			}

			if($isoform == 1) {
				my @seq_type_locus = ($columns[3], $columns[4]);
				$seq_type{$gene_name}{$columns[2]} = \@seq_type_locus;
			}

			### Update gene locus
			if($columns[3] < $gene_ann{$gene_name}[0]) {
				$gene_ann{$gene_name}[0] = $columns[3];
				$seq_type{$gene_name}{"gene"}[0]= $columns[3];
			}
			if($columns[4] > $gene_ann{$gene_name}[1]) {
				$gene_ann{$gene_name}[1] = $columns[4];
				$seq_type{$gene_name}{"gene"}[1]= $columns[4];
			}
		}
	}
}


### Load chromosome
print STDERR "Reading reference sequence file\n";
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
print STDERR "Get coding sequence\n";
foreach my $gene_name ( sort keys %coding_ann ) {
	$gene_ann{$gene_name}[3] = substr($chr_seq, $gene_ann{$gene_name}[0] - 1, $gene_ann{$gene_name}[1] - $gene_ann{$gene_name}[0] + 1);

	foreach my $start ( sort {$a <=> $b} keys %{$coding_ann{$gene_name}} ) {
		$coding_ann{$gene_name}{$start}[3] = substr($chr_seq, $start - 1, $coding_ann{$gene_name}{$start}[1] - $start + 1);
	}
}

### Add gene SNP and Indel annotation
print STDERR "Calculating annotation\n";
foreach my $gene_name ( sort keys %gene_ann ) {

	print STDERR "Processing gene: $gene_name\n";

	# Coding SNP protein change
	if( exists $coding_ann{$gene_name} ) {
		my $gene = new GeneSNPlist($snps);
		
		my $results = $gene->get_gene_snps( $chromosome, $gene_ann{$gene_name}[0], $gene_ann{$gene_name}[1],
						$gene_ann{$gene_name}[2], $gene_name, 1, $gene_ann{$gene_name}[3], %{$coding_ann{$gene_name}});
		
		$gene->get_protein_changes();
	}

	# Calculate SNP and Indel annotation type
	foreach my $seq_type ( keys %{$seq_type{$gene_name}} ) {
		my $seq_type_start = $seq_type{$gene_name}{$seq_type}[0];
		my $seq_type_end   = $seq_type{$gene_name}{$seq_type}[1];

		# SNPs
		foreach my $snp_pos ( sort {$a <=> $b} keys %{$snps->{snps}} ) {

			if( ($snp_pos >= $seq_type_start) && ($snp_pos <= $seq_type_end) ) {

				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$snps->{snps}{$snp_pos}{stype} = $seq_type;
				}

				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $snps->{snps}{$snp_pos}{stype} ne "CDS") {
						$snps->{snps}{$snp_pos}{stype} = $seq_type;
						$snps->{snps}{$snp_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(	($snps->{snps}{$snp_pos}{stype} ne "CDS") && 
						($snps->{snps}{$snp_pos}{stype} ne "exon")
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
				( ($del_pos >= $seq_type_start) && ($del_pos <= $seq_type_end) ) #||
				#( ($deletions->{indels}{$del_pos}{end} >= $seq_type_start) && ($deletions->{indels}{$del_pos}{end} <= $seq_type_end) )
			) {
					
				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$deletions->{indels}{$del_pos}{stype} = "CDS";
					$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
				}

				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $deletions->{indels}{$del_pos}{stype} ne "CDS") {
						$deletions->{indels}{$del_pos}{stype} = $seq_type;
						$deletions->{indels}{$del_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(     ($deletions->{indels}{$del_pos}{stype} ne "CDS") &&
						($deletions->{indels}{$del_pos}{stype} ne "exon")
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
				( ($ins_pos >= $seq_type_start) && ($ins_pos <= $seq_type_end) ) #||
				#( ($insertions->{indels}{$ins_pos}{end} >= $seq_type_start) && ($insertions->{indels}{$ins_pos}{end} <= $seq_type_end) )
			) {
					
				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$insertions->{indels}{$ins_pos}{stype} = $seq_type;
					$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
				}


				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $insertions->{indels}{$ins_pos}{stype} ne "CDS") {
						$insertions->{indels}{$ins_pos}{stype} = $seq_type;
						$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(	($insertions->{indels}{$ins_pos}{stype} ne "CDS") &&
						($insertions->{indels}{$ins_pos}{stype} ne "exon")
					){
						$insertions->{indels}{$ins_pos}{stype} = "intronic";
						$insertions->{indels}{$ins_pos}{gene_id} = $gene_name;
					}
				}
			}
		}

		# HDR
		foreach my $hdr_beg ( sort {$a <=> $b} keys %{$HDRs->{hdrs}} ) {
			my $hdr_end = $HDRs->{hdrs}{$hdr_beg}{end};

			for(my $i = $hdr_beg; $i <= $hdr_end; $i++) {
					
				if( ($i >= $seq_type_start) && ($i <= $seq_type_end) ) {
						
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$HDRs->{hdrs}{$hdr_beg}{stype} = "CDS";
						$HDRs->{hdrs}{$hdr_beg}{cds_absence}{$gene_name} = 1;
					}

					# Set sequence type to UTR
					elsif( $seq_type eq "exon" ) {
						if( $HDRs->{hdrs}{$hdr_beg}{stype} ne "CDS" ) {
							$HDRs->{hdrs}{$hdr_beg}{stype} = $seq_type;
						}
						$HDRs->{hdrs}{$hdr_beg}{utr_absence}{$gene_name} = 1;
					}

					# Set sequence type to intronic
					elsif( $seq_type eq "gene" ) {
						if(	($HDRs->{hdrs}{$hdr_beg}{stype} ne "CDS") &&
							($HDRs->{hdrs}{$hdr_beg}{stype} ne "exon")
						){
							$HDRs->{hdrs}{$hdr_beg}{stype} = "intronic";
						}
						$HDRs->{hdrs}{$hdr_beg}{gene_absence}{$gene_name} = 1;
					}
				}
			}
		}
	}

	print STDERR "Finished processing gene: $gene_name\n";
}

print STDERR "Finished calculating annotation, start printing results\n";


### Print SNP results
open SNPOUT, ">$outfolder/snp.annotation.Chr$chromosome.txt" or die "Cannot open SNP output file\n";

foreach my $pos (sort {$snps->{snps}{$a}{position} <=> $snps->{snps}{$b}{position}} keys %{$snps->{snps}} ) {
	
	print SNPOUT "$chromosome\t$pos\t" . $snps->{snps}{$pos}{ref_base} ."\t". $snps->{snps}{$pos}{new_base} ."\t". 
			$snps->{snps}{$pos}{quality} ."\t". $snps->{snps}{$pos}{support} ."\t". $snps->{snps}{$pos}{concordance} ."\t". 
			$snps->{snps}{$pos}{repetitive} ."\t". $snps->{snps}{$pos}{stype};

	if($snps->{snps}{$pos}{gene_id} ne "NA") {
		print SNPOUT "\t" . $snps->{snps}{$pos}{gene_id} . "\t1";
	}

	if($snps->{snps}{$pos}{cds_pos} != 0) {
		my $syn_nonsyn = "Syn";
		if($snps->{snps}{$pos}{ns_change} == 1) { $syn_nonsyn = "Nonsyn"; }

		print SNPOUT "\t" . $snps->{snps}{$pos}{cds_pos} . "\t" . $snps->{snps}{$pos}{codon_pos} . "\t$syn_nonsyn\t" . 
			$snps->{snps}{$pos}{ref_aa} . "\t" . $snps->{snps}{$pos}{new_aa};
	}
	print SNPOUT "\n";
}
close SNPOUT;


### Print deletion results
if($del_file ne "") {

	open DELOUT, ">$outfolder/del.annotation.Chr$chromosome.txt" or die "Cannot open deletion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$deletions->{indels}} ) {
		
		print DELOUT "$chromosome\t$start\t" . $deletions->{indels}{$start}{end} ."\t". $deletions->{indels}{$start}{seq} ."\t". 
				$deletions->{indels}{$start}{support} ."\t". $deletions->{indels}{$start}{concordance} ."\t".
				$deletions->{indels}{$start}{repetitive} ."\t". $deletions->{indels}{$start}{stype};

		if( $deletions->{indels}{$start}{gene_id} ne "NA") {
			print DELOUT "\t". $deletions->{indels}{$start}{gene_id};
		}
	
		print DELOUT "\n";
	}

	close DELOUT;
}


### Print insertion results
if($ins_file ne "") {
	open INSOUT, ">$outfolder/ins.annotation.Chr$chromosome.txt" or die "Cannot open insertion output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$insertions->{indels}} ) {

		print INSOUT "$chromosome\t$start\t" . $insertions->{indels}{$start}{end} ."\t". $insertions->{indels}{$start}{seq} ."\t". 
				$insertions->{indels}{$start}{support} ."\t". $insertions->{indels}{$start}{concordance} ."\t".
				$insertions->{indels}{$start}{repetitive} ."\t". $insertions->{indels}{$start}{stype};

		if( $insertions->{indels}{$start}{gene_id} ne "NA") {
			print INSOUT "\t". $insertions->{indels}{$start}{gene_id};
		}

		print INSOUT "\n";
	}
	
	close INSOUT;
}


### Print HDR results
if($hdr_file ne "") {
	open HDROUT, ">$outfolder/hdr.annotation.Chr$chromosome.txt" or die "Cannot open HDR output file\n";

	foreach my $start (sort {$a <=> $b} keys %{$HDRs->{hdrs}} ) {

		print HDROUT "$chromosome\t$start\t" . $HDRs->{hdrs}{$start}{end} ."\t". $HDRs->{indels}{$start}{dtype} ."\t". $HDRs->{hdrs}{$start}{stype};

		my $gene_list = "";
		my $utr_list = "";
		my $cds_list = "";

		foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$start}{gene_absence}} ) {
			$gene_list .= $gene_id . ",";
		}
		foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$start}{utr_absence}} ) {
			$utr_list .= $gene_id . ",";
		}
		foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$start}{cds_absence}} ) {
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

	close HDROUT;
}

print STDERR "Variant annotation finished\n";

exit(0);


### Read command line parameters
sub GetCom {
  my @usage = ("\nUsage: $0

Variation files and format (mandatory):
--format   STRING     Supported formats: SAM or shore
--snp      STRING     SNP file


Variation files in Shore format (optional):
--del      STRING     Deletion file
--ins      STRING     Insertion file
--hdr      STRING     HDR file (e.g. SVs, CNVs or unseq)


Functional SNP and indel annotation (mandatory):
--out      STRING     Output folder for result files
--chrom    STRING     Chromosome
--genome   STRING     Reference sequence file (chromosome names have to be equal to SNP file)
--gff      STRING     Gene annotation in GFF format

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "format=s", "snp=s", "del=s", "ins=s", "hdr=s", "out=s", "chrom=s", "genome=s", "gff=s");

	die("Please specify variant input format\n") unless defined($CMD{format});
	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify outfolder\n") unless defined($CMD{out});
	die("Please specify reference genome file\n") unless defined($CMD{genome});
	die("Please specify annotation gff file\n") unless defined($CMD{gff});
	die("Please specify chromosome of target region\n") unless defined($CMD{chrom});

	$format      = $CMD{format};
	$snp_file    = $CMD{snp};
	$outfolder   = $CMD{out};
	$chromosome  = $CMD{chrom};
	$refseq_file = $CMD{genome};
	$gff         = $CMD{gff};

	if(defined $CMD{del}) { $del_file = $CMD{del}; }
	if(defined $CMD{ins}) { $del_file = $CMD{ins}; }
	if(defined $CMD{hdr}) { $del_file = $CMD{hdr}; }
}

