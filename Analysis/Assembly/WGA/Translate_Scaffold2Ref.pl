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
#  Module: Analysis::Assembly::WGA::Translate_Scaffold2Ref.pl
#  Purpose:
#  In:
#  Out:
#




### User params
my $usage = "$0 ecotype assembly reference transin chromosome maxDivBases outfolder\n";

my $ecotype   = shift or die $usage;
my $assembly  = shift or die $usage;
my $reference = shift or die $usage;
my $transin   = shift or die $usage;
my $chr       = shift or die $usage;
my $maxDivBp  = shift or die $usage;
my $outfolder = shift or die $usage;


### Init other variables
my $chr_int = substr($chr, 3);

### Create folders
mkdir "$outfolder" or die;
mkdir "$outfolder/Alignments" or die;

### Run NUCMER
system("nucmer --mum -b $maxDivBp -g 90 -l 35 -c 80 -f --prefix=$outfolder/wgaMum $reference $assembly");

### Run delta filter
system("delta-filter -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");


### Parse reference file
open REFERENCE, $reference or die "Cannot open $reference\n";

my $id = "";
my $REF = "";

while(<REFERENCE>) {
        chomp;
	
	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		if($id eq $chr) {
			$REF .= $_;
		}
	}
}


### Parse assembly file
open ASSEMBLY, $assembly or die "Cannot open $assembly\n";

$id = "";
my %CTG = ();

while(<ASSEMBLY>) {
	chomp;

	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		$CTG{$id} .= $_;
	}
}
close ASSEMBLY;


### Parse translation input file
open TRANSIN, $transin or die "Cannot open $transin\n";

my %trans = ();
while(<TRANSIN>) {
	chomp;
	my @e = split ("\t", $_);
	if(exists $CTG{$e[0]}) {
		$trans{$e[0] ."#". $e[1]} = $_;
		#print "D1:$e[0]\t$e[1]\n";
	}
}

### Create alignments for each scaffold ID
open OUT, ">$outfolder/$chr.translation.txt" or die "Cannot open $outfolder/$chr.translation.txt\n";
foreach $id (keys %CTG) {

	### Create and open alignment for one scaffold
	system("show-aligns -r $outfolder/wgaMum.filter $chr $id > $outfolder/Alignments/alignments_$id.txt");
	open ALIGNMENT, "$outfolder/Alignments/alignments_$id.txt" or die "Cannot open $outfolder/Alignments/alignments_$id.txt\n";


	### Init variables
	my $state = 0;
	my $junk_line = "";
	my $HSPcounter = 0;
	my %HSP = ();


	### Parse alignment file of one scaffold
	while(<ALIGNMENT>) {
		chomp;

		### New HSP found ----------------------------------------------------------------------
		if( $_ =~ /BEGIN/ ) {
			$HSPcounter++;
			my @e = split(" ", $_);
			my %current_HSP = ( 'refbeg' => $e[5], 'refend' => $e[7], 'scaffbeg' => $e[10], 'scaffend' => $e[12] );
			$HSP{$HSPcounter} = \%current_HSP;

			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1;
		}


		### HSP finished, process information --------------------------------------------------
		elsif(  $_ =~ /END/ ) {

			### Check if translation input is overlapping
			foreach my $trans_locus (keys %trans) {
				my ($trans_id, $trans_pos) = split("#", $trans_locus);

				### Translate
				if( $id eq $trans_id ) {
					if( ($HSP{$HSPcounter}{scaffbeg} <= $trans_pos) && ($trans_pos <= $HSP{$HSPcounter}{scaffend}) ) {
						my $relative_distance = $trans_pos - $HSP{$HSPcounter}{scaffbeg};
						my $ref_pos = $HSP{$HSPcounter}{refbeg} + $relative_distance;
						print OUT "$ecotype\t$chr_int\t$ref_pos\t" . $trans{$trans_locus} . "\n";
					}
				}
			}

			### Reset
			$state = 0;
		}


		### Correct for stupid extra empty line after HSP-END ----------------------------------
		elsif( $_ eq "" ) {
			# Do nothing
		}


		### Reference line of alignment --------------------------------------------------------
		elsif( $state == 1 ) {
			$state = 2;
		}


		###Scaffold line of alignment ----------------------------------------------------------
		elsif( $state == 2 ) { 
			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1; 
		}
	}

	close ALIGNMENT;
}
close OUT;

print STDERR "Finished Scaffold Translation\n";
exit(0);
