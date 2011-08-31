#!/usr/bin/perl

# corrects indel alignments in map.list format
# indels are aligned most downstream
# improves overlap of different alignment tools (e.g. bwa or genomemapper)
# 

use strict;
use warnings;

my $usage = "\n$0 mapfile genomefile\n\n";
my $mapfile  = shift or die $usage;
my $genome   = shift or die $usage;


# Load genome fasta file
open GENOME, $genome or die "Cannot open input file\n";
my $chr = "";
my %gseq = ();
while( <GENOME> ) {
	chomp;
	if( $_ =~ />/ ) {
		$chr = substr($_, 1);
	}
	else {
		$gseq{$chr} .= $_;
	}
}


# Correct indels in map.list file
open MAPFILE, $mapfile or die "Cannot open input file\n";
while ( <MAPFILE> ) {
	my @e = split(/\t/, $_);

	# Check if alignment contains indels
	if( ($e[2] =~ /-/) || ($e[2] =~ /L/) ) {

		# Find indel start and end
		my $ali = $e[2];

		my $indelbeg = 0;
		my $indelend = 0;
		
		my $is_insertion = 0;
		my $is_deletion  = 0;
		
		my $split_read  = 0;
		my $multi_indel = 0;
		
		my $i = 0;
		
		# Aligned indels
		if( ($e[2] =~ /-/) && ($e[2] !~ /L/) ) {
			while( $i < length($ali) ) {

				# first indel position detected
				if( substr($ali, $i, 1) eq "-" ) {
				
					# Deletion
					if( substr($ali, $i-2, 1) eq "[" ) {
						$is_deletion = 1;
						$indelbeg = $i-2;
						$indelend = $i+1;
						$i += 2;
					}
					# Insertion
					elsif(substr($ali, $i-1, 1) eq "[" ) {
						$is_insertion = 1;
						$indelbeg = $i-1;
						$indelend = $i+2;
						$i += 3;
					}

					# Search end of indel
					while( $i < length($ali) ) {

						# Found consecutive variant
						if( substr($ali, $i, 1) eq "[" ) {
							
							# Variant is a SNP: Indel has ended
							if( (substr($ali, $i+1, 1) ne "-") && (substr($ali, $i+2, 1) ne "-") ) {
								$i += 4;
								last;
							}
							# Consecutive deletion
							elsif( $is_deletion && (substr($ali, $i+2, 1) eq "-") ) {
								$indelend += 4;
								$i += 4;
							}
							# Consecutive insertion
							elsif( $is_insertion && (substr($ali, $i+1, 1) eq "-") ) {
								$indelend += 4;
								$i += 4;
							}
							# complexe mix of insertion and deletion
							else {
								$multi_indel = 1;
								last;
							}
						}
						# Indel has ended
						else {
							last;
						}
					}

					# Check that no second indel is found in the read
					for(my $j = $i + 1; $j < length($ali); $j++ ) {
						if( substr($ali, $j, 1) eq "-" ) {
							$multi_indel = 1;
							last;
						}
					}

					last;
				}

				$i++;
			}
		}

		# Split reads
		elsif( $e[2] =~ /\[L\d+\]/ ) {
			$is_deletion = 1;
			$split_read  = 1;

			# Get begin and end of indel relative to read
			my $up = $`;
			$indelbeg = length($up);
			$indelend = $indelbeg + length($&) - 1;

			# Calculate real length of upstream sequence (not counting brackets/SNPs)
			my $uplen = 0;
			for(my $z = 0; $z < $indelbeg; $z++) {
				if(substr($up, $z, 1) eq "[") {
					$z+=3;
				}
				$uplen++;
			}

			# Get indel sequence from reference genome
			my $indellen = substr( $&, 2, length($&) - 3 );
			my $indelseq = substr( $gseq{$e[0]}, $e[1] + $uplen - 1, $indellen );

			# Create full indel alignment
			my $new_indelseq = "";
			for( my $x = 0; $x < length($indelseq); $x++) {
				$new_indelseq .= "[" . substr($indelseq, $x, 1) . "-]";
			}

			$ali = $up . $new_indelseq . substr($e[2], $indelend+1);
			$indelend = $indelbeg + length($new_indelseq) - 1;

			if( substr($ali, $i, 1) eq "-" ) {
				$multi_indel = 1;
			}

			#print $e[2] . "\t$indelbeg\t$indelend\t$indellen\t$e[1]\t$uplen\n";
			#print "$indelseq\t$new_indelseq\n";
			#print "$ali\t$indelbeg\t$indelend\n";
			#exit(0);
		}

		# Complex indel
		else {
			$multi_indel = 1;
		}



		#print "$indelbeg\t$indelend\t$multi_indel\t$is_deletion\t$is_insertion\t$ali\n";

		# Correct alignment
		if($multi_indel == 0) {
			my $upali    = substr($ali, 0, $indelbeg);
			my $indali   = substr($ali, $indelbeg, $indelend - $indelbeg + 1);
			my $downali  = substr($ali, $indelend+1); 
	
			my $notfixed = 1;
			my $n = $indelbeg;
	
			# Correct deletions
			if($is_deletion) {

				while( $notfixed && $indelbeg < length($ali) && $indelend < length($ali) ) {

					# First base of deletion is equal to first base after deletion: swap required
					if( substr($ali, $indelbeg+1, 1) eq substr($ali, $indelend+1, 1) ) {
					
						# First base in deletion is attached to the end of the upstream sequence (un-deleted)
						$upali .= substr($ali, $indelbeg+1, 1);
	
						# Downstream sequence is shortened by the first base, which is now part of the deletion
						$downali = substr($downali, 1);

						# Deletion swap: move deletion from fron to end
						$indali = substr($indali, 4) . substr($indali, 0, 4);

						# Update position of indel
						$indelbeg++;
						$indelend++;
						$ali = $upali . $indali . $downali;
					}
					# First base of deletion is NOT equal to first base after deletion: end correction procedure
					else {
						$notfixed = 0;
						last;
					}
				}
			
				# Print corrected alignment
				if($split_read == 1) {
					my $indellen = ($indelend - $indelbeg + 1) / 4;
					if($indellen > 30) {
						$ali = $upali . "[L" . $indellen . "]" . $downali;
					}
				}
				
				print $e[0]."\t".$e[1]."\t".$ali."\t".$e[3]."\t".$e[4]."\t".$e[5]."\t".$e[6]."\t".$e[7]."\t".$e[8]."\t".$e[9]."\t".$e[10];
			
			}
			# Correct insertions
			elsif($is_insertion) {
				# TODO
				print $_;
			}
			else {
				print STDERR "ERROR\n"; exit(0);
			}
		}
		else {
			print $_;
		}
	}
	# no indel alignment: keep original entry
	else {
		print $_;
	}
}

close MAPFILE;
close GENOME;

exit(0);
