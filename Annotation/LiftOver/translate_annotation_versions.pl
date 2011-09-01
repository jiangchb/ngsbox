#! /usr/bin/perl

use strict;

# translates positions in GFF files for different reference genome versions.
# optimized for bacterial genomes (e.g. Thermotoga)

my $usage = "$0  transtype<N1_LI | LI_N1 | N1_SB | SB_N1 | LI_SB | SB_LI>  infile_type<gff | snp | indel>  in_file  conversion_file\n";

my $transtype = shift or die $usage;	# N1_LI, LI_N1, N1_SB, SB_N1, LI_SB, SB_LI
my $file_type = shift or die $usage;	# gff, snp, indel
my $in_file   = shift or die $usage;
my $co_file   = shift or die $usage;


### Load conversion file
open CONV, $co_file or die "Cannot open $co_file\n";

my %trans = ();

while( <CONV> ) {
	chomp;
	my @a = split " ", $_;
	if   ( $transtype eq "N1_LI" ) {
		$trans{$a[0]} = $a[2];
	}
	elsif( $transtype eq "LI_N1" ) {
		$trans{$a[2]} = $a[0];
	}
	elsif( $transtype eq "N1_SB" ) {
		$trans{$a[0]} = $a[4];
	}
	elsif( $transtype eq "SB_N1" ) {
		$trans{$a[4]} = $a[0];
	}
	elsif( $transtype eq "LI_SB" ) {
		$trans{$a[2]} = $a[4];
	}
	elsif( $transtype eq "SB_LI" ) {
		$trans{$a[4]} = $a[2];
	}
}
close SNP;

open IN, $in_file or die "Cannot open $in_file\n";

while ( <IN> ) {
	chomp;
	
	### Ignore comments
	if( substr($_, 0, 1) eq "#" ) {
		print "$_\n";
	}

	### Translate
	else {

		### Translate GFF
		if( $file_type eq "gff") {
			my @a = split(/\t/, $_);
			print $a[0] ."\t". $a[1] ."\t". $a[2] ."\t". $trans{$a[3]} ."\t". $trans{$a[4]} ."\t". $a[5] ."\t". $a[6] ."\t". $a[7] ."\t". $a[8] ."\n";
		}

		### Translate homozygous_snp file
		elsif( $file_type eq "snp") {
			my @a = split(/\t/, $_);
			print $a[0] ."\t". $a[1] ."\t". $trans{$a[2]} ."\t". $a[3] ."\t". $a[4] ."\t". $a[5] ."\t". $a[6] ."\t". $a[7] ."\t". $a[8] ."\t". $a[9] ."\n";
		}

		### Translate indel file
		elsif( $file_type eq "indel") {
			my @a = split(/\t/, $_);
			print $a[0] ."\t". $a[1] ."\t". $trans{$a[2]} ."\t". $trans{$a[3]} ."\t". $a[4] ."\t". $a[5] ."\t". $a[6] ."\t". $a[7] ."\t". $a[8] ."\t". $a[9] ."\n";
		}
	}
}

close IN;
exit(0);
