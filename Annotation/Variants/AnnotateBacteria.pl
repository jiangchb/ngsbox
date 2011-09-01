#!/usr/bin/perl -w

# --------------------------------------------------------------------
# NGSBox: Variant Annotation
# 
# Annotate SNPs, indels, SVs and CNVs against any eukaryotic gene
# annotation in GFF format
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
my $format      = "";
my $snp_file    = "";
my $del_file    = "";
my $ins_file    = "";
my $hdr_file    = "";
my $outfolder   = "";

print STDERR "Reading command line parameters\n";

### Get command line options
my %CMD;
GetCom();


### Get SNP lists
print STDERR "Reading SNP file\n";
my $snps = new SNPlist();
$snps->get($snp_file, $format);


### Get deletion list
print STDERR "Reading deletion file\n";
my $deletions = new IndelList();
if($del_file ne "") { $deletions->get($del_file, $format); }


### Get insertion list
print STDERR "Reading insertion file\n";
my $insertions = new IndelList();
if($ins_file ne "") { $insertions->get($ins_file, $format); }


### Get HDR list
print STDERR "Reading HDR file\n";
my $HDRs = new HDRList();
if($hdr_file ne "") { $HDRs->get($hdr_file, $format); }



### Functional analysis of SNPs and indels
my %coding_ann = ();
my %gene_ann = ();
my %seq_type = ();


print STDERR "Reading GFF file\n";

### Get gene annotation from gff
open GFF, $gff or die "Cannot open gff file\n";
while( <GFF> ) {
	chomp;

	if(substr($_, 0, 1) ne "#") {

		# GFF format: chr, source, seq_type, start, end, score, orientation, frame, description
		my @c = split("\t", $_);
		my $gene_name = $c[8];

		# Add new gene
		if( ! exists $gene_ann{$gene_name} ) {

			# New gene
			my @gene_locus = ($c[0], $c[3], $c[4], $c[6], "");
			$gene_ann{$gene_name} = \@gene_locus;

			# Sequence type (gene, CDS, exone ...)
			my %type = ();
			$seq_type{$gene_name} = \%type;

			# Sequence type: gene
			my @seq_type_locus = ($c[0], $c[3], $c[4]);
			$seq_type{$gene_name}{"gene"} = \@seq_type_locus;

			# Sequence type: CDS
			if(! exists $coding_ann{$gene_name}) {
				my %gene = ();
				$coding_ann{$gene_name} = \%gene;
			}

			my @cds = ($c[0], $c[3], $c[4], $c[6], "");
			$coding_ann{$gene_name}{$c[3]} = \@cds;

			my @seq_type_locus2 = ($c[0], $c[3], $c[4]);
			$seq_type{$gene_name}{"CDS"} = \@seq_type_locus2;
				
			# Sequence type: exon 
			my @seq_type_locus3 = ($c[0], $c[3], $c[4]);
			$seq_type{$gene_name}{"exon"} = \@seq_type_locus3;


			### Update gene locus
			if($c[3] < $gene_ann{$gene_name}[1]) {
				$gene_ann{$gene_name}[1] = $c[3];
				$seq_type{$gene_name}{"gene"}[1]= $c[3];
			}
			if($c[4] > $gene_ann{$gene_name}[2]) {
				$gene_ann{$gene_name}[2] = $c[4];
				$seq_type{$gene_name}{"gene"}[2]= $c[4];
			}
		}
	}
}

### Load chromosome
print STDERR "Reading reference sequence file\n";
my %ref_seq = ();
my $current_chr = "";
open GENOME, $refseq_file or die "Cannot open reference sequence file\n";
while( <GENOME> ) {
	chomp;

	if(substr($_, 0, 1)  eq ">") {
		my @a = split(" ", $_);
		$current_chr = substr($a[0], 1);
	}
	else {
		$ref_seq{$current_chr} .= $_;
	}
}


### Extract gene seq
print STDERR "Get coding sequence\n";
foreach my $gene_name ( sort keys %coding_ann ) {
	my $chr = $gene_ann{$gene_name}[0];

	$gene_ann{$gene_name}[4] = substr($ref_seq{$chr}, $gene_ann{$gene_name}[1] - 1, $gene_ann{$gene_name}[2] - $gene_ann{$gene_name}[1] + 1);

	foreach my $start ( sort {$a <=> $b} keys %{$coding_ann{$gene_name}} ) {
		$coding_ann{$gene_name}{$start}[4] = substr($ref_seq{$chr}, $start - 1, $coding_ann{$gene_name}{$start}[2] - $start + 1);
	}
}


### Add gene SNP and Indel annotation
print STDERR "Calculating annotation\n";
foreach my $gene_name ( sort keys %gene_ann ) {

	print STDERR "Processing gene: $gene_name\n";
	my $chr = $gene_ann{$gene_name}[0];

	# Coding SNP protein change
	if( exists $coding_ann{$gene_name} ) {

		my $gene = new GeneSNPlist($snps);

		my $results = $gene->get_gene_snps( $chr, $gene_ann{$gene_name}[1], $gene_ann{$gene_name}[2],
			$gene_ann{$gene_name}[3], $gene_name, 1, $gene_ann{$gene_name}[4], %{$coding_ann{$gene_name}});

		$gene->get_protein_changes();
	}

	# Calculate SNP and Indel annotation type
	foreach my $seq_type ( keys %{$seq_type{$gene_name}} ) {
		my $seq_type_start = $seq_type{$gene_name}{$seq_type}[1];
		my $seq_type_end   = $seq_type{$gene_name}{$seq_type}[2];

		# SNPs
		foreach my $snp_pos ( sort {$a <=> $b} keys %{$snps->{snps}{$chr}} ) {

			if( ($snp_pos >= $seq_type_start) && ($snp_pos <= $seq_type_end) ) {

				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$snps->{snps}{$chr}{$snp_pos}{stype} = $seq_type;
				}

				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $snps->{snps}{$chr}{$snp_pos}{stype} ne "CDS") {
						$snps->{snps}{$chr}{$snp_pos}{stype} = $seq_type;
						$snps->{snps}{$chr}{$snp_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(	($snps->{snps}{$chr}{$snp_pos}{stype} ne "CDS") && 
						($snps->{snps}{$chr}{$snp_pos}{stype} ne "exon")
					){
						$snps->{snps}{$chr}{$snp_pos}{stype} = "intronic";
						$snps->{snps}{$chr}{$snp_pos}{gene_id} = $gene_name;
					}
				}
			}
		}

		# Deletions
		foreach my $del_pos ( sort {$a <=> $b} keys %{$deletions->{indels}{$chr}} ) {
				
			if( ($del_pos >= $seq_type_start) && ($del_pos <= $seq_type_end) ) {
					
				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$deletions->{indels}{$chr}{$del_pos}{stype} = "CDS";
					$deletions->{indels}{$chr}{$del_pos}{gene_id} = $gene_name;

					# Indel position in CDS: works only for single exon genes as found in bacteria!
					if($gene_ann{$gene_name}[3] eq "+") {
						$deletions->{indels}{$chr}{$del_pos}{cds_pos} = $del_pos - $seq_type_start + 1;
					}
					elsif($gene_ann{$gene_name}[3] eq "-") {
						$deletions->{indels}{$chr}{$del_pos}{cds_pos} = $seq_type_end - $del_pos + 1;
					}
				}

				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $deletions->{indels}{$chr}{$del_pos}{stype} ne "CDS") {
						$deletions->{indels}{$chr}{$del_pos}{stype} = $seq_type;
						$deletions->{indels}{$chr}{$del_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(     ($deletions->{indels}{$chr}{$del_pos}{stype} ne "CDS") &&
						($deletions->{indels}{$chr}{$del_pos}{stype} ne "exon")
					){
						$deletions->{indels}{$chr}{$del_pos}{stype} = "intronic";
						$deletions->{indels}{$chr}{$del_pos}{gene_id} = $gene_name;
					}
				}
			}
		}

		# Insertions
		foreach my $ins_pos ( sort {$a <=> $b} keys %{$insertions->{indels}{$chr}} ) {
				
			if( ($ins_pos >= $seq_type_start) && ($ins_pos <= $seq_type_end) ) {
					
				# Set sequence type to CDS
				if( $seq_type eq "CDS" ) {
					$insertions->{indels}{$chr}{$ins_pos}{stype} = $seq_type;
					$insertions->{indels}{$chr}{$ins_pos}{gene_id} = $gene_name;

					# Indel position in CDS: works only for single exon genes as found in bacteria!
					if($gene_ann{$gene_name}[3] eq "+") {
						$insertions->{indels}{$chr}{$ins_pos}{cds_pos} = $ins_pos - $seq_type_start + 1;
					}
					elsif($gene_ann{$gene_name}[3] eq "-") {
						$insertions->{indels}{$chr}{$ins_pos}{cds_pos} = $seq_type_end - $ins_pos + 1;
					}
				}


				# Set sequence type to UTR
				elsif( $seq_type eq "exon" ) {
					if( $insertions->{indels}{$chr}{$ins_pos}{stype} ne "CDS") {
						$insertions->{indels}{$chr}{$ins_pos}{stype} = $seq_type;
						$insertions->{indels}{$chr}{$ins_pos}{gene_id} = $gene_name;
					}
				}

				# Set sequence type to intronic
				elsif( $seq_type eq "gene" ) {
					if(	($insertions->{indels}{$chr}{$ins_pos}{stype} ne "CDS") &&
						($insertions->{indels}{$chr}{$ins_pos}{stype} ne "exon")
					){
						$insertions->{indels}{$chr}{$ins_pos}{stype} = "intronic";
						$insertions->{indels}{$chr}{$ins_pos}{gene_id} = $gene_name;
					}
				}
			}
		}

		# HDR
		foreach my $hdr_beg ( sort {$a <=> $b} keys %{$HDRs->{hdrs}{$chr}} ) {
			my $hdr_end = $HDRs->{hdrs}{$chr}{$hdr_beg}{end};

			for(my $i = $hdr_beg; $i <= $hdr_end; $i++) {
					
				if( ($i >= $seq_type_start) && ($i <= $seq_type_end) ) {
						
					# Set sequence type to CDS
					if( $seq_type eq "CDS" ) {
						$HDRs->{hdrs}{$chr}{$hdr_beg}{stype} = "CDS";
						$HDRs->{hdrs}{$chr}{$hdr_beg}{cds_absence}{$gene_name} = 1;
					}

					# Set sequence type to UTR
					elsif( $seq_type eq "exon" ) {
						if( $HDRs->{hdrs}{$chr}{$hdr_beg}{stype} ne "CDS" ) {
							$HDRs->{hdrs}{$chr}{$hdr_beg}{stype} = $seq_type;
						}
						$HDRs->{hdrs}{$chr}{$hdr_beg}{utr_absence}{$gene_name} = 1;
					}

					# Set sequence type to intronic
					elsif( $seq_type eq "gene" ) {
						if(	($HDRs->{hdrs}{$chr}{$hdr_beg}{stype} ne "CDS") &&
							($HDRs->{hdrs}{$chr}{$hdr_beg}{stype} ne "exon")
						){
							$HDRs->{hdrs}{$chr}{$hdr_beg}{stype} = "intronic";
						}
						$HDRs->{hdrs}{$chr}{$hdr_beg}{gene_absence}{$gene_name} = 1;
					}
				}
			}
		}
	}

	print STDERR "Finished processing gene: $gene_name\n";
}

print STDERR "Finished calculating annotation, start printing results\n";


### Print results
foreach my $chr ( sort keys %ref_seq ) {

	### Print SNP results
	open SNPOUT, ">$outfolder/snp.annotation.txt" or die "Cannot open SNP output file\n";

	foreach my $pos (sort {$snps->{snps}{$chr}{$a}{position} <=> $snps->{snps}{$chr}{$b}{position}} keys %{$snps->{snps}{$chr}} ) {
	
		print SNPOUT "$chr\t$pos\t" . $snps->{snps}{$chr}{$pos}{ref_base} ."\t". $snps->{snps}{$chr}{$pos}{new_base} ."\t". 
			$snps->{snps}{$chr}{$pos}{quality} ."\t". $snps->{snps}{$chr}{$pos}{support} ."\t". $snps->{snps}{$chr}{$pos}{concordance} ."\t". 
			$snps->{snps}{$chr}{$pos}{repetitive} ."\t". $snps->{snps}{$chr}{$pos}{stype};

		if($snps->{snps}{$chr}{$pos}{gene_id} ne "NA") {
			my $gene_name = $snps->{snps}{$chr}{$pos}{gene_id};
			my $gene_length = $gene_ann{$gene_name}[2] - $gene_ann{$gene_name}[1] + 1;
			print SNPOUT "\t" . $snps->{snps}{$chr}{$pos}{gene_id} ."\t". $gene_length;
		}

		if($snps->{snps}{$chr}{$pos}{cds_pos} != 0) {
			my $syn_nonsyn = "Syn";
			if($snps->{snps}{$chr}{$pos}{ns_change} == 1) { $syn_nonsyn = "Nonsyn"; }

			print SNPOUT "\t" . $snps->{snps}{$chr}{$pos}{cds_pos} . "\t" . $snps->{snps}{$chr}{$pos}{codon_pos} . "\t$syn_nonsyn\t" . 
				$snps->{snps}{$chr}{$pos}{ref_aa} . "\t" . $snps->{snps}{$chr}{$pos}{new_aa};
		}
		print SNPOUT "\n";
	}
	close SNPOUT;


	### Print deletion results
	if($del_file ne "") {
	
		open DELOUT, ">$outfolder/del.annotation.txt" or die "Cannot open deletion output file\n";
		
		foreach my $start (sort {$a <=> $b} keys %{$deletions->{indels}{$chr}} ) {
		
			print DELOUT "$chr\t$start\t" . $deletions->{indels}{$chr}{$start}{end} ."\t". $deletions->{indels}{$chr}{$start}{seq} ."\t". 
				$deletions->{indels}{$chr}{$start}{support} ."\t". $deletions->{indels}{$chr}{$start}{concordance} ."\t".
				$deletions->{indels}{$chr}{$start}{repetitive} ."\t". $deletions->{indels}{$chr}{$start}{stype};

			if( $deletions->{indels}{$chr}{$start}{gene_id} ne "NA") {
				my $gene_name = $deletions->{indels}{$chr}{$start}{gene_id};
				my $gene_length = $gene_ann{$gene_name}[2] - $gene_ann{$gene_name}[1] + 1;
				print DELOUT "\t". $deletions->{indels}{$chr}{$start}{gene_id} ."\t". $gene_length;
			}

			if($deletions->{indels}{$chr}{$start}{cds_pos} != 0) {
				print DELOUT "\t". $deletions->{indels}{$chr}{$start}{cds_pos};
			}

			print DELOUT "\n";
		}

		close DELOUT;
	}


	### Print insertion results
	if($ins_file ne "") {
		open INSOUT, ">$outfolder/ins.annotation.txt" or die "Cannot open insertion output file\n";

		foreach my $start (sort {$a <=> $b} keys %{$insertions->{indels}{$chr}} ) {

			print INSOUT "$chr\t$start\t" . $insertions->{indels}{$chr}{$start}{end} ."\t". $insertions->{indels}{$chr}{$start}{seq} ."\t". 
				$insertions->{indels}{$chr}{$start}{support} ."\t". $insertions->{indels}{$chr}{$start}{concordance} ."\t".
				$insertions->{indels}{$chr}{$start}{repetitive} ."\t". $insertions->{indels}{$chr}{$start}{stype};

			if( $insertions->{indels}{$chr}{$start}{gene_id} ne "NA") {
				my $gene_name = $insertions->{indels}{$chr}{$start}{gene_id};
				my $gene_length = $gene_ann{$gene_name}[2] - $gene_ann{$gene_name}[1] + 1;
				print INSOUT "\t". $insertions->{indels}{$chr}{$start}{gene_id} ."\t". $gene_length;
			}

			if($insertions->{indels}{$chr}{$start}{cds_pos} != 0) {
				print INSOUT "\t". $insertions->{indels}{$chr}{$start}{cds_pos};
			}

			print INSOUT "\n";
		}
	
		close INSOUT;
	}


	### Print HDR results
	if($hdr_file ne "") {
		open HDROUT, ">$outfolder/hdr.annotation.txt" or die "Cannot open HDR output file\n";

		foreach my $start (sort {$a <=> $b} keys %{$HDRs->{hdrs}{$chr}} ) {

			print HDROUT "$chr\t$start\t" . $HDRs->{hdrs}{$chr}{$start}{end} ."\t". $HDRs->{hdrs}{$chr}{$start}{stype};

			my $gene_list = "";
			my $utr_list = "";
			my $cds_list = "";

			foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$chr}{$start}{gene_absence}} ) {
				$gene_list .= $gene_id . ",";
			}
			foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$chr}{$start}{utr_absence}} ) {
				$utr_list .= $gene_id . ",";
			}
			foreach my $gene_id ( sort keys %{$HDRs->{hdrs}{$chr}{$start}{cds_absence}} ) {
				$cds_list .= $gene_id . ",";
			}

			if($gene_list eq "") { $gene_list = "NA"; }
			else { 
				chop($gene_list); 
				$gene_list =~ s/ /_/g;
			}

			if($utr_list eq "") { $utr_list = "NA"; }
			else { 
				chop($utr_list); 
				$utr_list =~ s/ /_/g;
			}

			if($cds_list eq "") { $cds_list = "NA"; }
			else { 
				chop($cds_list);
				$cds_list =~ s/ /_/g;
			}

			#print HDROUT "\t$gene_list\t$utr_list\t$cds_list";
			print HDROUT "\t$cds_list";

			print HDROUT "\n";
		}

		close HDROUT;
	}
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
--genome   STRING     Reference sequence file (chromosome names have to be equal to SNP file)
--gff      STRING     Gene annotation in GFF format

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "format=s", "snp=s", "del=s", "ins=s", "hdr=s", "out=s", "genome=s", "gff=s");

	die("Please specify variant input format\n") unless defined($CMD{format});
	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify outfolder\n") unless defined($CMD{out});
	die("Please specify reference genome file\n") unless defined($CMD{genome});
	die("Please specify annotation gff file\n") unless defined($CMD{gff});

	$format      = $CMD{format};
	$snp_file    = $CMD{snp};
	$outfolder   = $CMD{out};
	$refseq_file = $CMD{genome};
	$gff         = $CMD{gff};

	if(defined $CMD{del}) { $del_file = $CMD{del}; }
	if(defined $CMD{ins}) { $ins_file = $CMD{ins}; }
	if(defined $CMD{hdr}) { $hdr_file = $CMD{hdr}; }
}

