#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use Indel;

package IndelList;

### Constructor
sub new {
	my $self = {
		format       => "",
		indel_file   => "",
		indels       => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $indel_file, $format) = @_;
	$self->{format}       = $format;
	$self->{indel_file}   = $indel_file;

	open INDELFILE, $indel_file or die "Cannot open indel file: " . $indel_file . "\n";

	### Shore indel format:  sample | chromosome | begin | end | length | seq | coreness | support | concordance | repetitiveness
	if($format eq "shore") {
		while (<INDELFILE>) {
			chomp;
			my ($sample, $chr, $begin, $end, $length, $seq, $coreness, $support, $concordance, $repetitiveness) = split "\t";

			$self->{indels}{$chr}{$begin} = new Indel;
			$self->{indels}{$chr}{$begin}->init( $sample, $chr, $begin, $end, $seq, $support, $concordance, $repetitiveness);
		}
	}

	### VCF format: chromosome | begin | refseq | indelseq | consensus_phred | indel_quality | mapping_quality | support | called_bases | base_quality
	elsif($format eq "vcf") {
		while (<INDELFILE>) {
			chomp;

			my ($chr, $begin,  $refseq, $indelseq, $consensus_phred, $indel_quality, $mapping_quality, $support, $allele_1, $allele_2, $junk1, $junk2, $junk3, $junk4, $junk5) = split "\t";


			# Transform indel sequence
			my @e = split("/", $indelseq);
			my $end = $begin + 1;
			my $seq = "";
			my $type = "";

			if($e[0] eq "*") {
				if(substr($e[1], 0, 1) eq "-") {
					$seq = substr($e[1], 1);
					$begin++;
					$end = $begin + length($seq) - 1;
					$type = "DEL";
				}
				elsif(substr($e[1], 0, 1) eq "+") {
					$seq = substr($e[1], 1);
					$type = "INS";
				}
			}
			else {
				if(substr($e[0], 0, 1) eq "-") {
					$seq = substr($e[0], 1);
					$begin++;
					$end = $begin + length($seq) - 1;
					$type = "DEL";
				}
				elsif(substr($e[0], 0, 1) eq "+") {
					$seq = substr($e[0], 1);
					$type = "INS";
				}
			}

			# Store
			$self->{indels}{$chr}{$begin} = new Indel;
			$self->{indels}{$chr}{$begin}->init( $type, $chr, $begin, $end, $seq, $support, $indel_quality, $mapping_quality );
		}
	}


	close INDELFILE;

	return(1);
}

1;
