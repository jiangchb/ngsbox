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
#  Module: Parser::FASTQ::trim_length.pl
#  Purpose: shorten reads in fastq files
#  In: 
#  Out:
#


my $usage  = "$0 minlength file1 file2\n";
my $min    = shift or die $usage;
my $file1   = shift or die $usage;
my $file2   = shift or die $usage;

my $adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC";

my @spacer_fwd = ("TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG", "TCGTATAACT", "TCTCGTATGC");
my @spacer_rev = ("CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA", "CGTAATAACT", "AGTTATACGA");

my %spacer = ();
$spacer{fwd} = \@spacer_fwd;
$spacer{rev} = \@spacer_rev;

open F1, $file1 or die "Cannot open input read 1 file\n";
open F2, $file2 or die "Cannot open input read 2 file\n";

open O1, ">$file1.filtered" or die "Cannot open output read 1 file\n";
open O2, ">$file2.filtered" or die "Cannot open output read 2 file\n";


while( <F1> ) {
        my $sh1   = $_;
        my $seq1  = <F1>;
        my $qh1   = <F1>;
        my $qual1 = <F1>;

        my $sh2   = <F2>;
        my $seq2  = <F2>;
        my $qh2   = <F2>;
        my $qual2 = <F2>;

	chomp $seq1;
	chomp $seq2;
	chomp $qual1;
	chomp $qual2;

	my $spacer_found = 0;

	# Illumina adapter found
	if( ($seq1 =~ m/$adapter/) || ($seq2 =~ m/$adapter/) ) {
		$spacer_found = 1;
	}


	foreach my $strand ( sort keys %spacer ) {

		my $spacer     = $spacer{$strand}[0];
		my $spacer_beg = $spacer{$strand}[1];
		my $spacer_end = $spacer{$strand}[2];


		### Read 1, cases 1/2/3:

		# start of spacer sequence identified
		if( (! $spacer_found) && ($seq1 =~ m/$spacer_beg/) ) {
		
			$spacer_found = 1;

			# Not case 1:
			if( length($`) >= $min) {
				my $print_seq1  = substr($seq1, 0, length($`));
				my $print_qual1 = substr($qual1, 0, length($`));
	
				print O1 "$sh1$print_seq1\n$qh1$print_qual1\n";
				print O2 "$sh2$seq2\n$qh2$qual2\n";
				#print "Case: $` $& $'\n";
			}
		}

		# end of spacer sequence identified
		if( (! $spacer_found) && ($seq1 =~ m/$spacer_end/) ) {
	
			$spacer_found = 1;
	
			my $remaining_length = length($`) + length($&) - length($spacer);
			
			# Not case 1:
			if( $remaining_length >= $min) {
				my $print_seq1  = substr($seq1, 0, $remaining_length);
				my $print_qual1 = substr($qual1, 0, $remaining_length);
	
				print O1 "$sh1$print_seq1\n$qh1$print_qual1\n";
				print O2 "$sh2$seq2\n$qh2$qual2\n";
				#print "Case: $` $& $'\n";
			}
		}

		### Read 2, cases 5/6/7:
		# start of spacer sequence identified
		if( (! $spacer_found) && ($seq2 =~ m/($spacer_beg)/) ) {

			$spacer_found = 1;

			# Not case 1:
			if( length($`) >= $min) {
				my $print_seq2  = substr($seq2, 0, length($`));
				my $print_qual2 = substr($qual2, 0, length($`));
	
				print O1 "$sh1$seq1\n$qh1$qual1\n";
				print O2 "$sh2$print_seq2\n$qh2$print_qual2\n";
				#print "Case: $` $& $'\n";
			}
		}

		# end of spacer sequence identified
		if( (! $spacer_found) && ($seq2 =~ m/($spacer_end)/) ) {
	
			$spacer_found = 1;
	
			my $remaining_length = length($`) + length($&) - length($spacer);
	
			# Not case 1:
			if( $remaining_length >= $min) {
				my $print_seq2  = substr($seq2, 0, $remaining_length);
				my $print_qual2 = substr($qual2, 0, $remaining_length);
	
				print O1 "$sh1$seq1\n$qh1$qual1\n";
				print O2 "$sh2$print_seq2\n$qh2$print_qual2\n";
				#print "Case: $` $& $'\n";
			}

		}
	}

	# Case 4 (good)
	if(! $spacer_found) {

		print O1 "$sh1$seq1\n$qh1$qual1\n";
		print O2 "$sh2$seq2\n$qh2$qual2\n";
	
	}
}

close F1;
close F2;
close O1;
close O2;

sub compare_seq {

	my ($seq1, $seq2) = @_;

	my $len = min(length($seq1), length($seq2));

	my $match = 0;
	my $mismatch = 0;

	for(my $i = 0; $i < $len; $i++) {
		if( 
			(substr($seq1, $i, 1) ne substr($seq2, $i, 1)) &&
			(substr($seq1, $i, 1) ne "N") &&
			(substr($seq2, $i, 1) ne "N")
		 ) {
			$mismatch++;
		}
	}

	if( $mismatch/$match > 0.1) {
		return(0);
	}
	else {
		return(1);
	}
}

sub min ($$) { $_[$_[0] > $_[1]] }


sub quality_filter {
	my ($seq, $qual) = @_;

}

exit(0);

### Jumping libraries and spacers
# Read: ====>
# Spacer: ---|******|--
# Cases:
#    |===============>                       <===============
# 1. |***|------------------------------------------------->|
# 2. |<----|*****|----------------------------------------->|
# 3. |<----------|*****|----------------------------------->|
# 4. |<------------------|******|-------------------------->|
# 5. |<------------------------------------|******|-------->|
# 6. |<------------------------------------------|******|-->|
# 7. |<-------------------------------------------------|***|
#
#
#print "Case1: $` $& $'\n";
