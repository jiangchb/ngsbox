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
#  Module: Analysis::Assembly::ConsensusAssembly::consensus_concatenation.pl
#  Purpose:
#  In:
#  Out:
#




### User parameters
my $usage = "$0 quality_reference_file quality_variant_file outfile centromeresfile\n";
my $r_file = shift or die $usage;
my $v_file = shift or die $usage;
my $o_file = shift or die $usage;
my $c_file = shift or die $usage;


### Files
open RFILE, $r_file or die "Cannot open file ".$r_file."\n";
open VFILE, $v_file or die "Cannot open file ".$v_file."\n";
open OFILE, ">".$o_file or die "Cannot open file ".$o_file."\n";
open CFILE, $c_file or die "Cannot open file ".$o_file."\n";


### Adjustable parameter
my $min_qual    = 15;
my $min_freq    = 0.7;
my $min_support = 2;
my $min_len     = 10;


### Statistics
my $contigs = 0;
my $max_contig = 0;


### Hashes
my %CTGlen = ();
my %QVAR   = ();
my %CENTRO = ();


### Read in quality var
while (my $l = <VFILE>) {
        my @a = split("\t", $l);
	if( ($a[7] >= $min_freq) and ($a[6] >= $min_support) and ($a[5] >= $min_qual) ) {
        	$QVAR{$a[1]."#".$a[2]} = $l;
	}
}
close VFILE;

### Read in centromeres
while (my $l = <CFILE>) {
        my @a = split("\t", $l);
	for (my $i = $a[1]; $i<= $a[2]; $i++) {
		$CENTRO{$a[0]."#".$i} = 1;
	}
}
close CFILE;

### Read in quality ref
my $last_chr = -1;
my $last_pos = -1;
my %CTG = ();

while (my $l = <RFILE>) {
	my @a = split "\t", $l;

	my $ref_chr = $a[1];
	my $ref_pos = $a[2];
	my $ref_nuc = $a[4];
	my $ref_q   = $a[5];

	if (not defined($CENTRO{$ref_chr."#".$ref_pos})) {

		# Chromosome change
		if($ref_chr != $last_chr) {
			$last_chr = $ref_chr;
			$last_pos = -1;
		}
	
	
		# Start new contig
		if( $last_pos == -1 ) {
			if( ($ref_q >= $min_qual) && ($a[7] >= $min_freq) && ($a[6] >= $min_support) ) {
				$contigs++;
				$CTG{$contigs} = ">$contigs | $ref_chr:$ref_pos#$ref_nuc";
				$last_pos = $ref_pos;
			}
		}
	
		# Contiguous ref call
		elsif( ($ref_chr == $last_chr) && ($ref_pos == $last_pos + 1) && ($ref_q >= $min_qual) && ($a[7] >= $min_freq) && ($a[6] >= $min_support) ) {
			$CTG{$contigs} .= $a[4];
			$last_pos = $ref_pos;
		}
	
		# Check if variant calls can fill the gap
		else {
			my $fill_flag = 1;
	
			for(my $i = $last_pos + 1; $i <= $ref_pos; $i++) {
	
				# Contiguous variant call found
				if( exists $QVAR{"$ref_chr#$i"} ) {
					my @a = split("\t", $QVAR{"$ref_chr#$i"});
	
					if( $a[4] ne "-" ) {
						$CTG{$contigs} .= $a[4];
					}
				}
	
				# Contig break if gap is not already closed
				else {
					# Gap is already closed
					if( ($ref_chr == $last_chr) && ($ref_pos == $i) && ($ref_q >= $min_qual) && ($a[7] >= $min_freq) && ($a[6] >= $min_support) ) {
						$CTG{$contigs} .= $ref_nuc;
					}
					# Contig break
					else {
						# ref base is used as first base of contig
						if( ($ref_q >= $min_qual) && ($a[7] >= $min_freq) && ($a[6] >= $min_support) ) {
							$contigs++;
							$CTG{$contigs} = ">$contigs | $ref_chr:$ref_pos#$ref_nuc";
							last;
						}
						# ref base quality is too low, have to look for next high quality ref base
						else {
							$last_pos = -1;
							$fill_flag = 0;
							last;
						}
					}
				}
			}
	
			if( $fill_flag == 1) {
				$last_pos = $ref_pos;
			}
		}
	}
	else {
		$last_pos = -1;
	}
}
close RFILE;

$contigs = 0;
foreach my $contig_nr (sort {$a<=>$b} keys %CTG) {
	my ($header, $seq) = split("#", $CTG{$contig_nr});
	my $ctg_len = length($seq);

	if( $ctg_len >= $min_len ) {
		$contigs++;
		if($ctg_len > $max_contig) {
			$max_contig = $ctg_len;
		}
		print OFILE "$header\n$seq\n";
	}
}


#####################################
# Print stats

print STDERR "Contigs: $contigs\n";
print STDERR "Contig max: $max_contig\n";

exit(0);

