#!/usr/bin/perl 

# --------------------------------------------------------------------
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use SNP;

package SNPlist;

### Constructor
sub new {
	my $self = {
		chromosome   => "",
		snp_file     => "",
		snps         => {},
	};
	bless $self;
	return $self;
}
		

sub get {

	### Init
	my ($self, $chromosome, $snp_file) = @_;
	$self->{chromosome}   = $chromosome;
	$self->{snp_file}     = $snp_file;

	open SNPFILE, $snp_file or die "Cannot open snp file: " . $snp_file . "\n";


	### WGA SNP format:  sample | chromosome | position | scaffold | scaffold position | refbase | snp
	while (<SNPFILE>) {
		chomp;

		my ($sample, $chr, $position, $scaff, $scaff_pos, $refbase, $snp) = split "\t";

		if ( $chr == $chromosome ) {
			$self->{snps}{$position} = new SNP();
			$self->{snps}{$position}->init( $sample, $chromosome, $position, $refbase, $snp );
		}
	}
	close SNPFILE;

	return(1);
}

1;
