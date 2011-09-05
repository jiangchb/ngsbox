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
#  Module: Annotation::Variants::SNPlist.pm
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use SNP;

package SNPlist;

### Constructor
sub new {
	my $self = {
		format       => "",
		snp_file     => "",
		snps         => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $snp_file, $format) = @_;
	$self->{format}   = $format;
	$self->{snp_file} = $snp_file;

	open SNPFILE, $snp_file or die "Cannot open snp file: " . $snp_file . "\n";


	### SHORE quality variant format: sample | chromosome | position | refbase | snp | quality | support | concordance | repetitiveness
	if($format eq "shore") {
		while (<SNPFILE>) {
			chomp;
			my ($sample, $chr, $position, $refbase, $snp, $quality, $support, $concordance, $repetitiveness) = split "\t";

			$self->{snps}{$chr}{$position} = new SNP();
			$self->{snps}{$chr}{$position}->init( $sample, $chr, $position, $refbase, $snp, $quality, $support, $concordance, $repetitiveness );
		}
	}

	### VCF format: chromosome | position | refbase | snp | consensus_phred | SNP_quality | mapping_quality | support | called_bases | base_quality
	elsif($format eq "vcf") {

		my %IUB = (
			M => 'AC', 
			R => 'AG', 
			W => 'AT', 
			S => 'CG', 
			Y => 'CT', 
			K => 'GT', 
		);

		while (<SNPFILE>) {
			chomp;

			my ($chr, $position, $refbase, $iub, $consensus_phred, $SNP_quality, $mapping_quality, $support, $called_bases, $base_quality) = split "\t";
			
			$refbase = uc($refbase);
			$iub = uc($iub);


			# Transform SNP format
			my $snp = "";
				
			if( $iub !~ /[ACGTacgt]/ ) {

				if( substr($IUB{$iub}, 0, 1) ne $refbase ) {
					$snp = substr($IUB{$iub}, 0, 1);
				}
				else {
					$snp = substr($IUB{$iub}, 1, 1);
				}
			}
			else {
				$snp = $iub;
			}

			# Store
			$self->{snps}{$chr}{$position} = new SNP();
			$self->{snps}{$chr}{$position}->init( "NA", $chr, $position, $refbase, $snp, $SNP_quality, $support, $consensus_phred, $mapping_quality );
		}
	}

	close SNPFILE;

	return(1);
}

1;
