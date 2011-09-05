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
#  Module: Parser::SHORE::Variants::make_polymorphic_genome.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 reffasta qualityvariation\n\n";

my $fasta = shift or die $usage;
my $var = shift or die $usage;

open FILE, $fasta or die $usage;

my %SIZE = ();

my %REF = ();
my $id = "";
my $seq = "";

################################################
# Read in reference fasta file
################################################
while (my $line = <FILE>) {
	if (substr($line, 0, 1) eq ">" ) {
		if ($id ne "") {
			print STDERR "Got >$id<\n";
			for (my $i = 0; $i < length($seq); $i++) {
				my $pos = $i + 1;
				$REF{$id."#".$pos} = substr($seq, $i, 1);
			}
			$SIZE{$id} = length($seq);			
		}
		my @a = split " ", $line;
		$id = substr($a[0], 1, length($a[0])-1);
		$seq = "";
	}
	else {
		chomp($line);
		$seq .= $line;
	}
}
if ($id ne "") {
	print STDERR "Got >$id<\n";
	for (my $i = 0; $i < length($seq); $i++) {
        	my $pos = $i + 1;
                $REF{$id."#".$pos} = substr($seq, $i, 1);
        }
        $SIZE{$id} = length($seq);
}
close FILE;

foreach my $key (keys %SIZE) {
	print STDERR ">", $key, "<\t>", $SIZE{$key}, "<\n";
}

################################################
# Read in SNPs and indels
################################################
open FILE, $var or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	#if ($a[4] ne "-") {
	if ($a[5] > 24) {     #### NEED TO CHECK IF THIS IS REALLY THE QUALITY COLUMN!!!
		$REF{$a[1]."#".$a[2]} = $a[4];
	}
	#}
	#elsif ($a[4] eq "-") {
#		#for (my $i = 0; $i < $a[4]; $i++) {
#		#my $pos = $a[3] + $i;
#		$REF{$a[2]."#".$pos} = substr($REF{$a[2]."#".$pos}, 1, length($REF{$a[2]."#".$pos})-1);
#		#}
#	}
#	else {
#		my $pos = $a[3] - 1;
 #               $REF{$a[2]."#".$pos} .= $a[5];
#	}
}

################################################
# Output
################################################
open OUT, "> new_genome.fa";
foreach my $chr (sort {$a <=> $b} keys %SIZE) {
	print STDERR "Outputting $chr\n";
	print OUT ">$chr\n";
	my $num = 0;
	for (my $i = 1; $i <= $SIZE{$chr}; $i++) {
		for (my $k = 0; $k < length($REF{$chr."#".$i}); $k++) {
			my $c = substr($REF{$chr."#".$i}, $k, 1);
			if ($c ne "-") {
				$num++;
				print OUT $c;
				print OUT "\n" if $num%79 == 0;
			}
		}

		#$num++;
		#print OUT $REF{$chr."#".$i};
		#print OUT "\n" if $num%79 == 0;
	}
	print OUT "\n";
}



