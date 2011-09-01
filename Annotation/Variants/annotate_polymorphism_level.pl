#!/usr/bin/perl
use strict;
use warnings;

my $usage = "$0 AnnotationGFF quality_reference.txt quality_variants.txt\n";
my $anno = shift or die $usage;
my $reference = shift or die $usage;
my $variants = shift or die $usage;

my %REF = ();
my %POLY = ();
my %GENES = ();
my %ANNOTATION = ();
my %EXON = ();
my %UTR = ();

###############################################
## Read in files
get_anno_file(\%GENES, \%EXON, \%UTR, \%ANNOTATION, $anno);
print STDERR "Got anno file\n";
get_qual_file(\%REF, $reference);
print STDERR "Got reference file\n";
get_qual_file(\%POLY, $variants);
print STDERR "Got variants file\n";


GENE: foreach my $tair (sort keys %GENES) {

	#print STDERR ">$tair<",$GENES{$tair},"\n";

	my @a = split "#", $GENES{$tair};
	my $chr = $a[0];
	my $begin = $a[1];
	my $end = $a[2];

	###########################################
	# Check for ref call
	my $c = 0;

	my $c_cons = 0;
	my $c_poly = 0;

	my $c_exon = 0;
	my $c_exon_cons = 0;
	my $c_exon_poly = 0;

	my $c_intron = 0;
	my $c_intron_cons = 0;
	my $c_intron_poly = 0;

	my $c_utr = 0;
	my $c_utr_cons = 0;
	my $c_utr_poly = 0;

	for (my $i = $begin; $i <= $end; $i++) {
		my $cons = 0;
		my $poly = 0;
		$c++;
		# Conservered of polymorph
		if (defined($REF{$chr."#".$i})) {
			$c_cons++;
			$cons = 1;
		}
		elsif (defined($POLY{$chr."#".$i})) {
                        $c_poly++;
                        $poly = 1;
                }

		# By annotation
		if (defined($EXON{$chr."#".$i})) {
			$c_exon++;
                        if ($cons == 1) {
				$c_exon_cons++;
			}
			if ($poly == 1) {
                                $c_exon_poly++;
                        }
		}
		elsif (defined($UTR{$chr."#".$i})) {
                        $c_utr++;
                        if ($cons == 1) {
                                $c_utr_cons++;
                        }
			if ($poly == 1) {
                                $c_utr_poly++;
                        }
                }
		else  {
                        $c_intron++;
                        if ($cons == 1) {
                                $c_intron_cons++;
                        }
			if ($poly == 1) {
                                $c_intron_poly++;
                        }
                }
	}

	print $tair, "\t", $ANNOTATION{$tair}, "\t", $c, "\t", $c_cons, "\t", $c_poly, "\t";
	if (($c_cons+$c_poly) == 0) {
		print "-1";
	}
	else {
		printf("%.3f", $c_cons/($c_cons+$c_poly));
	}


	print "\t", $c_exon_cons, "\t", $c_exon_poly, "\t";
	if (($c_exon_cons+$c_exon_poly) == 0)  {
		print "-1";
	}
	else {
		printf("%.3f", $c_exon_cons/($c_exon_cons+$c_exon_poly));
	}

	print "\t", $c_intron_cons, "\t", $c_intron_poly, "\t";
	if (($c_intron_cons+$c_intron_poly) == 0) {
		print "-1";
	}
	else {
		printf("%.3f", $c_intron_cons/($c_intron_cons+$c_intron_poly));
	}

	print "\t", $c_utr_cons, "\t", $c_utr_poly, "\t";
	if (($c_utr_cons+$c_utr_poly) == 0) {
		print "-1";
	}
	else {
		printf("%.3f", $c_utr_cons/($c_utr_cons+$c_utr_poly));
	}	

	print "\n";

}






###################################################
# subroutines

sub get_qual_file {
	my ($ref, $file) = @_;
	
	open FILE, $file or die "Cannot open file: $file\n";
	while (my $line=<FILE>) {
		my @a = split " ", $line;
		if ($a[5] >= 25) {
			${$ref}{$a[1]."#".$a[2]} = $a[4]."#".$a[8];
		}
	}
	close FILE;
}

sub get_anno_file {
	my ($gene_ref, $exon_ref, $utr_ref, $anno_ref, $file) = @_;

	open FILE, $file or die "Cannot open file\n";

	while (my $line=<FILE>) {
		my @a = split " ", $line;
		$a[0] =~ s/Chr//g;
		if ($a[2] eq "gene" or $a[2] eq "pseudogene") {
			my @b = split ";", $a[8];
			$b[0] =~ s/ID=//g;
			${$gene_ref}{$b[0]} = $a[0]."#".$a[3]."#".$a[4];
			${$anno_ref}{$b[0]} = $a[2];
		}
		if ($a[2] eq "CDS" or $a[2] eq "pseudogenic_exon") {
			for (my $i = $a[3]; $i <= $a[4]; $i++) {
				${$exon_ref}{$a[0]."#".$i} = 1;
			}
		}
		if ($a[2] eq "three_prime_UTR" or $a[2] eq "five_prime_UTR") {
                        for (my $i = $a[3]; $i <= $a[4]; $i++) {
                                ${$utr_ref}{$a[0]."#".$i} = 1;
                        }
                }
	}

	close FILE or die "file won't close\n";
}












