#!/usr/bin/perl

use strict;
use warnings;

my $fasta_masked   = shift or die;
my $fasta_unmasked = shift or die;

my $seq_masked   = "";
my $seq_unmasked = "";

### Read masked fasta file
open FASTA, $fasta_masked or die "Cannot open masked fasta file\n";
while (<FASTA>) {
	chomp;
	if( $_ =~ /^>/ ) { 
		# Do nothing
	}
	else {
		$seq_masked .= $_;
	}
}
close FASTA or die "Cannot close fasta file\n";


### Read unmasked fasta file
open FASTA, $fasta_unmasked or die "Cannot open unmasked fasta file\n";
while (<FASTA>) {
	chomp;
	if( $_ =~ /^>/ ) { 
		# Do nothing
	}
	else {
		$seq_unmasked .= $_;
	}
}
close FASTA or die "Cannot close unmasked fasta file\n";


### Compare
print ">diff_masked_unmasked\n";
for(my $i = 0; $i < length($seq_masked); $i++) {
	my $unmasked_base = substr($seq_unmasked, $i, 1);
	if( substr($seq_masked, $i, 1) ne $unmasked_base ) {
		print $unmasked_base;
	}
}

exit(0);
