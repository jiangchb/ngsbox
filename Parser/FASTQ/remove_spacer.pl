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
#  Module: Parser::FASTQ::JumpLO.pl
#  Purpose: removes spacer sequences from jumping libraries
#  In: read1 and read2 sequences as fastq files
#  Out: cleaned fastq files
#

use List::Util qw[min max];

my $usage  = "$0 seedlength minlength file1 file2\n";
my $seed    = shift or die $usage;
my $min     = shift or die $usage;
my $file1   = shift or die $usage;
my $file2   = shift or die $usage;

my $adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC";

my $spacer_fwd_string = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";
my $spacer_rev_string = "CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA";

my @spacer_fwd = ( $spacer_fwd_string, substr($spacer_fwd_string, 0, $seed), substr($spacer_fwd_string, length($spacer_fwd_string) - $seed, $seed) );
my @spacer_rev = ( $spacer_rev_string, substr($spacer_rev_string, 0, $seed), substr($spacer_rev_string, length($spacer_rev_string) - $seed, $seed) );
#print "$spacer_fwd[0]\t$spacer_fwd[1]\t$spacer_fwd[2]\n";
#print "$spacer_rev[0]\t$spacer_rev[1]\t$spacer_rev[2]\n";



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

	my $determined = 0;

	# Check Illumina filter flag
	my @header1 = split(":", $sh1);
	my @header2 = split(":", $sh2);
	if( ($header1[7] eq "Y") || ($header2[7] eq "Y") ) {
		$determined = 1;
		#print "$sh1$seq1\n$qh1$qual1\n$sh2$seq2\n$qh2$qual2\n";
		#print "$qual1\t$qual2\n";
	}

	# Check Illumina adapter sequence
	if( ($seq1 =~ m/$adapter/) || ($seq2 =~ m/$adapter/) ) {
		$determined = 1;
		#print "$sh1$seq1\n$qh1$qual1\n$sh2$seq2\n$qh2$qual2\n";
		#print "$seq1\t$seq2\n";
	}


	foreach my $strand ( sort keys %spacer ) {

		my $spacer_seq = $spacer{$strand}[0];
		my $spacer_beg = $spacer{$strand}[1];
		my $spacer_end = $spacer{$strand}[2];


		### Read 1, cases 1/2/3:

		# start of spacer sequence identified
		if( (! $determined) && ($seq1 =~ m/$spacer_beg/) ) {
			my $spacer_in_read = $& . $';
		
			if( &compare_seq($spacer_seq, $spacer_in_read) ) {
				$determined = 1;

				# Not case 1:
				if( length($`) >= $min) {
					my $print_seq1  = substr($seq1, 0, length($`));
					my $print_qual1 = substr($qual1, 0, length($`));
		
					print O1 "$sh1$print_seq1\n$qh1$print_qual1\n";
					print O2 "$sh2$seq2\n$qh2$qual2\n";

					#print "$` $& $' " . length($`) . "\t" . length($spacer_seq) . "\t$spacer_seq\t$spacer_in_read\n";
				}
				#else {
				#	print "$` $& $' " . length($`) . "\t" . length($spacer_seq) . "\t$spacer_seq\t$spacer_in_read\n";
				#}
			}
		}

		# end of spacer sequence identified
		if( (! $determined) && ($seq1 =~ m/.*($spacer_end)/) ) {

			# define position and sequence of spacer in read
			my $remaining_read = max(0, length($&) - length($spacer_seq));
			my $spacer_in_read = substr($seq1, $remaining_read, min(length($&), length($spacer_seq)) );

			# define region of spacer for alignment
			my $cut_spacer = max(0, length($spacer_seq) - length($&));
			my $remaining_spacer = substr($spacer_seq, $cut_spacer);

			if( &compare_seq($remaining_spacer, $spacer_in_read) ) {
				$determined = 1;	

				# Not case 1:
				if( $remaining_read >= $min) {
					my $print_seq1  = substr($seq1, 0, $remaining_read);
					my $print_qual1 = substr($qual1, 0, $remaining_read);
		
					print O1 "$sh1$print_seq1\n$qh1$print_qual1\n";
					print O2 "$sh2$seq2\n$qh2$qual2\n";

					#print "$& $'\t" . length($&) . "\t$cut_spacer\t$remaining_read\t$remaining_spacer\t$spacer_in_read\n";
				}
				#else {
				#	print "$& $'\t" . length($&) . "\t$cut_spacer\t$remaining_read\t$remaining_spacer\t$spacer_in_read\n";
				#}
			}
			
		}

		### Read 2, cases 5/6/7:
		# start of spacer sequence identified
		if( (! $determined) && ($seq2 =~ m/($spacer_beg)/) ) {
			my $spacer_in_read = $& . $';

			if( &compare_seq($spacer_seq, $spacer_in_read) ) {
				$determined = 1;

				# Not case 1:
				if( length($`) >= $min) {
					my $print_seq2  = substr($seq2, 0, length($`));
					my $print_qual2 = substr($qual2, 0, length($`));
		
					print O1 "$sh1$seq1\n$qh1$qual1\n";
					print O2 "$sh2$print_seq2\n$qh2$print_qual2\n";

					#print "$` $& $' " . length($`) . "\t" . length($spacer_seq) . "\t$spacer_seq\t$spacer_in_read\n";
				}
				#else {
				#	print "$` $& $' " . length($`) . "\t" . length($spacer_seq) . "\t$spacer_seq\t$spacer_in_read\n";
				#}
			}
		}

		# end of spacer sequence identified
		if( (! $determined) && ($seq2 =~ m/.*($spacer_end)/) ) {

			# define position and sequence of spacer in read
			my $remaining_read = max(0, length($&) - length($spacer_seq));
			my $spacer_in_read = substr($seq2, $remaining_read, min(length($&), length($spacer_seq)) );

			# define region of spacer for alignment
			my $cut_spacer = max(0, length($spacer_seq) - length($&));
			my $remaining_spacer = substr($spacer_seq, $cut_spacer);

			if( &compare_seq($remaining_spacer, $spacer_in_read) ) {
				$determined = 1;
	
				# Not case 1:
				if( $remaining_read >= $min) {
					my $print_seq2  = substr($seq2, 0, $remaining_read);
					my $print_qual2 = substr($qual2, 0, $remaining_read);
	
					print O1 "$sh1$seq1\n$qh1$qual1\n";
					print O2 "$sh2$print_seq2\n$qh2$print_qual2\n";

					#print "$& $'\t" . length($&) . "\t$cut_spacer\t$remaining_read\t$remaining_spacer\t$spacer_in_read\n";
				}
				#else {
				#	print "$& $'\t" . length($&) . "\t$cut_spacer\t$remaining_read\t$remaining_spacer\t$spacer_in_read\n";
				#}
			}
		}
	}

	# Case 4 (good)
	if(! $determined) {

		print O1 "$sh1$seq1\n$qh1$qual1\n";
		print O2 "$sh2$seq2\n$qh2$qual2\n";
		#print "$sh1$sh2";
	}
}

close F1;
close F2;
close O1;
close O2;

sub compare_seq {

	my ($seq1, $seq2) = @_;

	my $len = min(length($seq1), length($seq2));

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

	my $ratio = $mismatch/$len;

	if( $ratio > 0.2) {
		return(0);
	}
	else {
		return(1);
	}
}

# sub min ($$) { $_[$_[0] > $_[1]] }


# sub max ($$) { $_[$_[0] < $_[1]] }


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
##if( (! $determined) && ($seq1 =~ m/$spacer_end(?!.*$spacer_end)/) ) {
